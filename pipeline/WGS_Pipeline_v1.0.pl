#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my $usage=<<'USAGE';
    
WGS data process pipeline of hg38 for Sentieon Platform (DIPSEQ T1).
Version 1 : 2021-10-28
Author : luoshuzhen, CANCER INSTITUTE
    
Options:
    -list            STR  reads file info list. <require>
    -outdir          STR  the output directory. <require>
    -split           STR  split commands. [default:true]
    -h                    Help Information

!!! Note : each sample must has only one sequencing library !!!!!
USAGE

my ($List,$outdir,$split);

GetOptions(
    'list=s' => \$List,
    'outdir=s' => \$outdir,
    'split=s' => \$split,
    'h|?' => sub{die "$usage\n";},
);

die "$usage\n" unless (defined($List) && defined($outdir));

$split ||= "true";
my $platform="COMPLETE_GENOMICS"; ## sequencing platform, fixed[COMPLETE_GENOMICS]
my $dir_log="$outdir/Run_WGS_Log";  ## recording samples
die "The directory for sample recording $dir_log does not exist!\n" if(!-e $dir_log);

## database 
my $REF_DIR = "/opt/bioinfo/database/gatk_bundle/hg38";
my $reference = "$REF_DIR/Homo_sapiens_assembly38.fasta";
my $Dbsnp = "$REF_DIR/dbsnp_151.hg38.vcf.gz";
my $kgIndel = "$REF_DIR/1000G_phase1.snps.high_confidence.hg38.vcf.gz";
my $millIndel = "$REF_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz";

## Parameters of SOAPnuke
my $SOAPnuke = "/opt/bioinfo/tools/SOAPnuke2.1.6/SOAPnuke"; ## Version: 2.1.6
my $SOAPnukePara = "-l 5 -q 0.5 -n 0.1 -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG -T 44";

## Sentieon
my $Sentieon = "/opt/bioinfo/tools/sentieon-genomics-202010.02/bin/sentieon";
my $SentieonCPU = 48;
my $Server = "export SENTIEON_LICENSE=sentieon:8990";

## tools
#my $perl = "/Share_Data/Opt/miniconda3/bin/perl";       ## perl-5.26.2
my $perl = "/usr/bin/perl";       ## v5.16.3
my $soapnuke_stat = "/opt/bioinfo/tools/bin/soapnuke2_stat.pl";
my $fqcheck = "/opt/bioinfo/tools/bin/fqcheck33";

## check for reference file
die "The reference $reference does not exist!\n" if(!-e $reference);
die "The reference $Dbsnp does not exist!\n" if(!-e $Dbsnp);
die "The reference $kgIndel does not exist!\n" if(!-e $kgIndel);
die "The reference $millIndel does not exist!\n" if(!-e $millIndel);

## file list for each sample
my %sampleFile;
my %seqLib;   ## sequencing library for each sample must = 1;
open LI,$List or die $!;
while(<LI>){
    chomp;
    my ($File,$Lib,$Sample,$Chip) = (split/\t/,$_)[1,3,4,6];
    if($File =~ /1.fq.gz$/){
	$File =~ s/\_1.fq.gz//;
    }elsif($File =~ /2.fq.gz$/){
	$File =~ s/\_2.fq.gz//;
    }
    $sampleFile{$Sample}{$Chip}{$File} = 1;
    $seqLib{$Sample} = $Lib;
}
close LI;

## fastq clean and alignment by sentiton
for my $samp (sort keys %sampleFile){
    
    `mkdir -p $outdir/$samp`;
    `mkdir -p $outdir/$samp/fq_raw`;
    `mkdir -p $outdir/$samp/clean`;
    `mkdir -p $outdir/$samp/sentieon`;
    `mkdir -p $outdir/$samp/shell`;
    
    `touch $outdir/$samp/fq_raw/fqstq_list.log`;  ## fqstq list for check 

    my @bamarray;       ## multiple bams
    
    ## merge and clean fastq
    my $cleanCmd="## merge fastq file, cat fastq file, ALT-Aware Mapping for each DIPSEQT1 id\n";
    $cleanCmd.="echo ========== merge and clean start at : `date` ========== && \\\n";

    ## alt-aware mapping
    my $alnCmd="echo ========== read alignment by sentieon start at : `date` ========== && \\\n";

    ## creat scrips for fqcheck
    my $checkCmd="echo ========== fqcheck start at : `date` ========== && \\\n";
    
    for my $chip (sort keys %{$sampleFile{$samp}}){
	`mkdir -p $outdir/$samp/clean/$chip`;

	my @fq1s;  
	my @fq2s;
	for my $file (sort keys %{$sampleFile{$samp}{$chip}}){
	    my $lane=(split/\//,$file)[-1];
	    my $fq1="$file\_1.fq.gz";
	    my $fq2="$file\_2.fq.gz";

	    ## merge fqstqs for each chip id
	    $cleanCmd.="if [ -e \"$fq1\" ];then echo -e \"$fq1\" >> $outdir/$samp/fq_raw/fqstq_list.log;fi && \\\n"; 
	    $cleanCmd.="if [ -e \"$fq2\" ];then echo -e \"$fq2\" >> $outdir/$samp/fq_raw/fqstq_list.log;fi && \\\n";
	    push @fq1s,$fq1;
	    push @fq2s,$fq2;
	}
	
	my $fq1_list=join(" ",@fq1s);
	my $fq2_list=join(" ",@fq2s);
	my $fq1_num=@fq1s;
	my $fq2_num=@fq2s;
	my $fqout1="$outdir/$samp/clean/$chip/$chip\_1.clean.fq.gz";    ## output fq1 name for SOAPnuke
	my $fqout2="$outdir/$samp/clean/$chip/$chip\_2.clean.fq.gz";    ## output fq2 name for SOAPnuke
	
	## cat fq for each DNBSEQ chip id
	if($fq1_num > 1 ){
	    $cleanCmd.="cat $fq1_list > $outdir/$samp/fq_raw/$chip\_1.fq.gz && \\\n";
	}else{
	    $cleanCmd.="ln -s $fq1_list $outdir/$samp/fq_raw/$chip\_1.fq.gz && \\\n";
	}
	if($fq2_num > 1 ){
	    $cleanCmd.="cat $fq2_list > $outdir/$samp/fq_raw/$chip\_2.fq.gz && \\\n";
	}else{
	   $cleanCmd.="ln -s $fq2_list $outdir/$samp/fq_raw/$chip\_2.fq.gz && \\\n";
	}
	
	## gzip test compressed file integrity
	$cleanCmd.="gzip -t $outdir/$samp/fq_raw/$chip\_1.fq.gz && \\\n";
	$cleanCmd.="gzip -t $outdir/$samp/fq_raw/$chip\_2.fq.gz && \\\n";
	
	## fastq Clean by SOAPnuke
	$cleanCmd.="$SOAPnuke filter -1 $outdir/$samp/fq_raw/$chip\_1.fq.gz -2 $outdir/$samp/fq_raw/$chip\_2.fq.gz $SOAPnukePara -o $outdir/$samp/clean/$chip -C $chip\_1.clean.fq.gz -D $chip\_2.clean.fq.gz && \\\n";
	$cleanCmd.="perl $soapnuke_stat $outdir/$samp/clean/$chip/Basic_Statistics_of_Sequencing_Quality.txt $outdir/$samp/clean/$chip/Statistics_of_Filtered_Reads.txt > $outdir/$samp/clean/$chip.clean.stat && \\\n";
	$cleanCmd.="if [ -e \"$fqout1\" ];then rm -rf $outdir/$samp/fq_raw/$chip\_1.fq.gz;fi && \\\n";
        $cleanCmd.="if [ -e \"$fqout2\" ];then rm -rf $outdir/$samp/fq_raw/$chip\_2.fq.gz;fi && \\\n";

	$checkCmd.="$fqcheck -r $fqout1 -c $outdir/$samp/clean/$chip\_1.clean.fqcheck && \\\n";
	$checkCmd.="$fqcheck -r $fqout2 -c $outdir/$samp/clean/$chip\_2.clean.fqcheck && \\\n";

	## clean reads align by Sentieon (BWA-mem)
	my $lib=$seqLib{$samp};
	$alnCmd.="$Sentieon bwa mem -M -t $SentieonCPU -R \"\@RG\\tID:$chip\\tLB:$lib\\tPU:$chip\\tSM:$samp\\tPL:$platform\\tCN:BGI\" $reference $fqout1 $fqout2 |$Sentieon util sort -r $reference -o $outdir/$samp/sentieon/$chip.sorted.bam -t $SentieonCPU --sam2bam -i - && \\\n";
	push @bamarray,"$outdir/$samp/sentieon/$chip.sorted.bam";
    }
    $cleanCmd.="echo ========== download and clean fastq end at : `date` ========== && \\\n";
    $alnCmd.="echo ========== read alignment by sentieon end at : `date` ========== && \\\n";
    
    ## Merge, indel-realignment, BQSR and QC by Sentieon
    ## merge multiple bams and mark PCR duplicates
    my $bqsrCmd="echo ========== merge and markdup start at : `date` ========== && \\\n";
    my $instring = join " -i ", @bamarray;
    $bqsrCmd.="$Sentieon driver -t $SentieonCPU -i $instring --algo LocusCollector --fun score_info $outdir/$samp/sentieon/$samp.mkdup.score.txt.gz && \\\n";
    $bqsrCmd.="$Sentieon driver -t $SentieonCPU -i $instring --algo Dedup --score_info $outdir/$samp/sentieon/$samp.mkdup.score.txt.gz --metrics $outdir/$samp/sentieon/$samp.mkdup.metrics.txt $outdir/$samp/sentieon/$samp.mkdup.bam && \\\n";
    $bqsrCmd.="if [ -e \"$outdir/$samp/sentieon/$samp.mkdup.bam\" ];then rm -rf $outdir/$samp/sentieon/*.sorted.bam*;fi && \\\n";
    $bqsrCmd.="echo ========== merge and mkdrdup end at : `date` ========== && \\\n";

    ## indel-realignment
    $bqsrCmd.="echo ========== indel-realignment start at : `date` ========== && \\\n";
    $bqsrCmd.="$Sentieon driver -t $SentieonCPU -r $reference -i $outdir/$samp/sentieon/$samp.mkdup.bam --algo Realigner -k $kgIndel -k $millIndel $outdir/$samp/sentieon/$samp.realn.bam && \\\n";
    $bqsrCmd.="echo ========== indel-realignment end at : `date` ========== && \\\n";

    ## bqsr & snp, indel calling
    $bqsrCmd.="echo ==========  BQSR start at : `date` ========== && \\\n";
    $bqsrCmd.="$Sentieon driver -t $SentieonCPU -r $reference -i $outdir/$samp/sentieon/$samp.realn.bam --algo QualCal -k $Dbsnp -k $kgIndel -k $millIndel $outdir/$samp/sentieon/$samp.bqsr.RECAL_PRE.txt && \\\n";
    $bqsrCmd.="$Sentieon driver -t $SentieonCPU -r $reference -i $outdir/$samp/sentieon/$samp.realn.bam -q $outdir/$samp/sentieon/$samp.bqsr.RECAL_PRE.txt --algo QualCal -k $Dbsnp -k $kgIndel -k $millIndel $outdir/$samp/sentieon/$samp.bqsr.RECAL_POST.txt --algo ReadWriter $outdir/$samp/sentieon/$samp.bqsr.bam --algo Haplotyper -d $Dbsnp $outdir/$samp/sentieon/$samp.Haplotyper.vcf.gz && \\\n";
    $bqsrCmd.="if [ -e \"$outdir/$samp/sentieon/$samp.bqsr.bam\" ];then rm -rf $outdir/$samp/sentieon/$samp.mkdup.bam $outdir/$samp/sentieon/$samp.mkdup.bam.bai $outdir/$samp/sentieon/$samp.realn.bam $outdir/$samp/sentieon/$samp.realn.bam.bai;fi && \\\n";
    $bqsrCmd.="echo ==========  BQSR end at : `date` ========== && \\\n";

    ## qc by Sentieon 
    my $qcCmd="echo ========== qc start at : `date` ========== && \\\n";
    $qcCmd.="$Sentieon driver -t $SentieonCPU -r $reference -i $outdir/$samp/sentieon/$samp.bqsr.bam --algo GCBias --summary $outdir/$samp/sentieon/$samp.bqsr.total.GC_SUMMARY.txt $outdir/$samp/sentieon/$samp.bqsr.total.GC_METRIC.txt --algo MeanQualityByCycle $outdir/$samp/sentieon/$samp.bqsr.total.MQ_METRIC.txt --algo QualDistribution $outdir/$samp/sentieon/$samp.bqsr.total.QD_METRIC.txt --algo InsertSizeMetricAlgo $outdir/$samp/sentieon/$samp.bqsr.total.IS_METRIC.txt --algo AlignmentStat $outdir/$samp/sentieon/$samp.bqsr.total.ALN_METRIC.txt --algo WgsMetricsAlgo $outdir/$samp/sentieon/$samp.bqsr.total.WGS_METRIC.txt --algo SequenceArtifactMetricsAlgo --dbsnp $Dbsnp $outdir/$samp/sentieon/$samp.bqsr.total.ARTIFACT_METRICS && \\\n";
    $qcCmd.="$Sentieon plot GCBias -o $outdir/$samp/sentieon/$samp.bqsr.total.GC_METRIC.pdf $outdir/$samp/sentieon/$samp.bqsr.total.GC_METRIC.txt && \\\n";
    $qcCmd.="$Sentieon plot MeanQualityByCycle -o $outdir/$samp/sentieon/$samp.bqsr.total.MQ_METRIC.pdf $outdir/$samp/sentieon/$samp.bqsr.total.MQ_METRIC.txt && \\\n";
    $qcCmd.="$Sentieon plot QualDistribution -o $outdir/$samp/sentieon/$samp.bqsr.total.QD_METRIC.pdf $outdir/$samp/sentieon/$samp.bqsr.total.QD_METRIC.txt && \\\n";
    $qcCmd.="$Sentieon plot InsertSizeMetricAlgo -o $outdir/$samp/sentieon/$samp.bqsr.total.IS_METRIC.pdf $outdir/$samp/sentieon/$samp.bqsr.total.IS_METRIC.txt && \\\n";
    $qcCmd.="echo ========== qc end at : `date` ========== && \\\n";

    ## output commands
    if($split =~ /true/i){
        my $check="#!/bin/bash\nset -euo pipefail\n";
        my $cmdEnd="echo -e SUCCESS\n";
        &Print("$outdir/$samp/shell/s1_clean\_$samp.sh",$check.$cleanCmd.$cmdEnd);
        &Print("$outdir/$samp/shell/s2_aln\_$samp.sh",$check.$Server."\n".$alnCmd.$cmdEnd);
        &Print("$outdir/$samp/shell/s3_bqsr\_$samp.sh",$check.$Server."\n".$bqsrCmd.$cmdEnd);
        &Print("$outdir/$samp/shell/s4_qc\_$samp.sh",$check.$Server."\n".$qcCmd.$cmdEnd);
        $checkCmd.="echo ========== fqcheck end at : `date` ========== && \\\n";
        &Print("$outdir/$samp/shell/s5_fqcheck\_$samp.sh",$check.$checkCmd.$cmdEnd);
    }else{
        my $finish_log="$dir_log/finish\_$samp.sign";
        open COM,">$outdir/$samp/shell/run\_$samp.sh" || die $!;
        print COM "#!/bin/bash\n";
        print COM "set -euo pipefail\n";
        print COM $cleanCmd,"echo -e \"clean fastq is OK\" > $finish_log && \\\n########################################\n";
        print COM $Server,"\n";
        print COM $alnCmd,"echo -e \"alignment is OK\" >> $finish_log && \\\n########################################\n";
        print COM $bqsrCmd,"echo -e \"merge and bqsr is OK\" >> $finish_log && \\\n########################################\n";
        print COM $qcCmd,"echo -e \"qc by sentieon is OK\" >> $finish_log && \\\n########################################\n";
        print COM "echo -e \"$samp is operation succeeded\" >> $finish_log && \\\n";
        close COM;
	open CK,">$outdir/$samp/shell/fqcheck\_$samp.sh" || die $!;
	print CK "#!/bin/bash\n";
	print CK $checkCmd;
	print CK "echo ==========  fqcheck end at : `date` ==========\n";
	close CK;
    }
}

sub Print{
    my($out,$cmd)=@_;
    open OUT,">$out" || die $!;
    print OUT $cmd;
    close OUT;
}

