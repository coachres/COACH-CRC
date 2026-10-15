#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd 'abs_path';
my $usage = <<'USAGE';
    
Analysis pipeline for WGBS data of DNBSEQ and illumina Platform.
Version 1.0 : 2023-01-11

Options:
	-list		STR  reads info list.  <require>
	-config		STR  configuration of parameters to run the pipeline. <require>
	-outdir		STR  the output directory. <require>
	-rmCleanReads	STR  whether remove clean reads after finish BitMapperBS mapping. default is true. [optional: true|false]
	-h		    Help Information
                            		 
!!! This pipeline is wrote base on chromosome name with "chr", be careful when replace default files!!!!!
USAGE

my($List,$config,$outdir,$rmCleanReads);

GetOptions(
	'list=s' => \$List,
	'config=s' => \$config,
	'outdir=s' => \$outdir,
	'rmCleanReads=s' => \$rmCleanReads,
	'h|?' => sub{die "$usage\n";},
);
die "$usage\n" unless (defined($List) && defined($config) && defined($outdir));
$outdir = abs_path($outdir);
$rmCleanReads ||= "true";
die "String input for '-rmCleanReads' is typo: must be true or false\n" unless ($rmCleanReads=~/true/i || $rmCleanReads=~/false/i);

my(%samples);

my($BIN, $DATAPATH, $reference, $ref_lambda, $MethylID, $BitMapperBS, $BitCPU, $MethylDackel, $MethylCPU, $MethylPara); ## basic
my($monitor, $queue, $priority, $subProjectName); ## monitor
my($length, $Use_SOAPnuke, $SOAPnuke, $SOAPnukePara, $mbias1, $mbias2, $fqcheck); # fq clean 
my($samtools, $sambamba, $sambambaCPU, $cgibed, $java, $picard, $Rbin); # dup, coverage, depth
my($Barcode_stat, $Sample_stat); # stat

################################ Read configuration file ###################################
open CFG,"$config" or die $!;
while(<CFG>){
	chomp;
	next if (/^#/);
	next unless (/\w+/);
	if($_ =~/^queue=(.*?);/){$queue=$1;}
	elsif($_ =~/^monitor=(.*?);/){$monitor=$1;}
	elsif($_ =~/^priority=(.*?);/){$priority=$1;}
	elsif($_ =~/^subProjectName=(.*?);/){$subProjectName=$1;}
	elsif($_ =~/^BitCPU=(.*?);/){$BitCPU=$1;}
	elsif($_ =~/^DATAPATH=(.*?);/){$DATAPATH=$1;}
	elsif($_ =~/^BitMapperBS=(.*?);/){$BitMapperBS=$1;}
	elsif($_ =~/^reference=(.*?);/){$reference=$1;}
	elsif($_ =~/^ref_lambda=(.*?);/){$ref_lambda=$1;}
	elsif($_ =~/^MethylID=(.*?);/){$MethylID=$1;}
	elsif($_ =~/^BIN=(.*?);/){$BIN=$1;}
	elsif($_ =~/^SOAPnuke=(.*?);/){$SOAPnuke=$1;}
	elsif($_ =~ /^cgibed=(.*?);/){$cgibed=$1;}
	elsif($_ =~/^samtools=(.*?);/){$samtools=$1;}
	elsif($_ =~/^sambamba=(.*?);/){$sambamba=$1;}
	elsif($_ =~/^sambambaCPU=(.*?);/){$sambambaCPU=$1;}
	elsif($_ =~/^java=(.*?);/){$java=$1;}
	elsif($_ =~/^picard=(.*?);/){$picard=$1;}
	elsif($_ =~/^Rbin=(.*?);/){$Rbin=$1;}
	elsif($_ =~/^MethylDackel=(.*?);/){$MethylDackel=$1;}
	elsif($_ =~/^MethylCPU=(.*?);/){$MethylCPU=$1;}
	elsif($_ =~/^Use_SOAPnuke=(.*?);/){$Use_SOAPnuke=$1;}
	elsif($_ =~/^length=(.*?);/){$length=$1;}
	elsif($_ =~/^mbias1=(.*?);/){$mbias1=$1;}
	elsif($_ =~/^mbias2=(.*?);/){$mbias2=$1;}
	elsif($_ =~ /^SOAPnukePara=(.*?);/){$SOAPnukePara=$1;}
}

################################### main parameters #####################################
$BIN ||= "/opt/bioinfo";
$DATAPATH ||= "/opt/bioinfo/database";
die "The bin directory $BIN does not exist!\n" unless(-e $BIN);
die "The database directory $DATAPATH does not exist!\n" unless(-e $DATAPATH);

## ref
$reference ||= "$DATAPATH/ucsc.hg38.fasta";
$ref_lambda ||= "$DATAPATH/lambda.fa";
$cgibed ||= "$DATAPATH/ucsc.hg38.cpgIslandExt.list";
$MethylID ||= "$DATAPATH/hg38_Methyl_site.ID.gz";
die "The reference $reference does not exist!\n" unless(-e $reference);
die "The reference index $reference.bai does not exist!\n" unless(-e "$reference.fai");
die "The ref_lambda $ref_lambda does not exist!\n" unless(-e $ref_lambda);
die "The ref_lambda index $ref_lambda.bai does not exist!\n" unless(-e "$ref_lambda.fai");
die "The cgi bed $cgibed does not exist!\n" unless(-e $cgibed);
die "The C site ID $MethylID does not exist!\n" unless(-e $MethylID);

## tools
$SOAPnuke ||= "$BIN/SOAPnuke";
$samtools ||= "$BIN/samtools";
$BitMapperBS ||= "$BIN/bitmapperBS";
$BitCPU ||= 8;
$sambamba ||= "$BIN/sambamba";
$sambambaCPU ||= 8;
$MethylDackel ||= "$BIN/MethylDackel";
$MethylCPU ||= 24;
$monitor ||= "$BIN/monitor";
$fqcheck ||= "$BIN/fqcheck33";
$Barcode_stat ||= "$BIN/Barcode_stat.pl";
$Sample_stat ||= "$BIN/Sample_stat.v2.pl";
$picard ||= "$BIN/picard.jar";
$Rbin ||= "/opt/R/4.1.1/bin";

## clean
unless(defined $length && $length ne ''){$length = 100}
unless(defined $mbias1 && $mbias1 ne ''){$mbias1 = 0}
unless(defined $mbias2 && $mbias2 ne ''){$mbias2 = 20}
$Use_SOAPnuke ||= "true";
$SOAPnukePara ||= "-l 5 -q 0.5 -n 0.1 -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG";

my @r = ('a'..'z','A'..'Z');
$subProjectName ||= join '', map { $r[int rand @r] } 0..6;

###################################Process sequence###################################
open LL,$List or die $!;
while(<LL>){
	chomp;
	my ($q,$id,$fq1,$fq2);
	my ($samp,$qual,$readsId,$lib,$lane,$lane_barcode,$fqprefix)=split /\s+/, $_;
	$samples{$samp}{$lib}{$lane_barcode} = 1;

	my $rgid = $lane_barcode;
	$rgid=~s/_read//;

	`mkdir -p $outdir/$samp/shell`;
	if($Use_SOAPnuke=~/true/i){
		`mkdir -p $outdir/$samp/01.Clean`;
	
		if($qual==33){$q=2;}
		else{$q=1;}
	
		if($readsId=~/old/i){$id=0;}
		elsif($readsId=~/new/i){$id=1;}
		else{die "Please offer reads id is old or new!";}

		my $rawfq1 = "";
		my $rawfq2 = "";
		if(-d $fqprefix){
			$rawfq1 = "$fqprefix/$lane_barcode\_1.fq.gz";
			$rawfq2 = "$fqprefix/$lane_barcode\_2.fq.gz";
		}
		else{
			$rawfq1 = "$fqprefix\_1.fq.gz";
			$rawfq2 = "$fqprefix\_2.fq.gz";
		}
		die "Do not exists $fqprefix*.fq.gz\n" unless (-e $rawfq1 && -e $rawfq2);
        
		$fq1 = "$outdir/$samp/01.Clean/$lane_barcode/$lane_barcode\_1.clean.fq.gz";
		$fq2 = "$outdir/$samp/01.Clean/$lane_barcode/$lane_barcode\_2.clean.fq.gz";
	
		open CLEAN,">$outdir/$samp/shell/Clean.$lane_barcode.sh" || die $!;
		print CLEAN "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` on `uname -a` ========== && \\\n";
		print CLEAN "$SOAPnuke filter -1 $rawfq1 -2 $rawfq2 $SOAPnukePara -t $mbias1,0,$mbias2,0 -Q $q -G $q -T 4 --seqType $id -o $outdir/$samp/01.Clean/$lane_barcode -C $lane_barcode\_1.clean.fq.gz -D $lane_barcode\_2.clean.fq.gz && \\\n";
		print CLEAN "/usr/bin/gzip -t $fq1 $fq2 && \\\n";
		print CLEAN "/usr/bin/perl $BIN/soapnuke2_stat.pl $outdir/$samp/01.Clean/$lane_barcode/Basic_Statistics_of_Sequencing_Quality.txt $outdir/$samp/01.Clean/$lane_barcode/Statistics_of_Filtered_Reads.txt > $outdir/$samp/01.Clean/$lane_barcode/$lane_barcode.clean.stat && \\\n";
		print CLEAN "echo ==========end at : `date` ========== && \\\n";
		print CLEAN "echo Still_waters_run_deep 1>&2\n";
		close CLEAN;
	}
	elsif($Use_SOAPnuke=~/false/i){
		$fq1 = "$fqprefix\_1.fq.gz";
		$fq2 = "$fqprefix\_2.fq.gz";
	}
	else{
		die "String input for 'Use_SOAPnuke' is typo: must be true or false";
	}

	open FASTQ,">$outdir/$samp/shell/fqcheck.$lane_barcode.sh" || die $!;
	print FASTQ "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` on `uname -a` ========== \n";
	print FASTQ "$fqcheck -r $fq1 -c $outdir/$samp/01.Clean/$lane_barcode/$lane_barcode\_1.clean.fqcheck && \\\n";
	print FASTQ "$fqcheck -r $fq2 -c $outdir/$samp/01.Clean/$lane_barcode/$lane_barcode\_2.clean.fqcheck && \\\n";
	print FASTQ "/usr/bin/perl $BIN/fqcheck_distribute.pl $outdir/$samp/01.Clean/$lane_barcode/$lane_barcode\_1.clean.fqcheck $outdir/$samp/01.Clean/$lane_barcode/$lane_barcode\_2.clean.fqcheck -o $outdir/$samp/01.Clean/$lane_barcode/$lane_barcode.clean. && \\\n";
	print FASTQ "rm -f $outdir/$samp/01.Clean/$lane_barcode/*.txt $outdir/$samp/01.Clean/$lane_barcode/log && \\\n";
	print FASTQ "echo ==========end at : `date` ========== && \\\n";
	print FASTQ "echo Still_waters_run_deep 1>&2\n";
	close FASTQ;
	
	`mkdir -p $outdir/$samp/02.BitBS`;
	open BITMAP,">$outdir/$samp/shell/Alignment.$lane_barcode.sh" || die $!;
	print BITMAP "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` on `uname -a` ========== && \\\n";
	print BITMAP "export LD_LIBRARY_PATH=/opt/bioinfo/BitMapperBS/v1.0.2.3-centOS7/SSE4.2/htslib_aim/lib && \\\n";
    print BITMAP "# $BitMapperBS --search $reference --seq1 $fq1 --seq2 $fq2 -t $BitCPU --bam --mapstats $outdir/$samp/02.BitBS/$lane_barcode.BitMapperBS.PE.stats | $samtools addreplacerg -@ $BitCPU -r '\@RG\\tID:$rgid\\tLB:$lib\\tPU:$lane\\tSM:$samp\\tCN:BGI' >$outdir/$samp/02.BitBS/$lane_barcode.bam - && \\\n";
    print BITMAP "mkdir -p /data/scratch/$lane_barcode && \\\n";
    print BITMAP "$samtools sort -@ 24 -m 10G --output-fmt BAM $outdir/$samp/02.BitBS/$lane_barcode.bam -o $outdir/$samp/02.BitBS/$lane_barcode.sort.bam -T /data/scratch/$lane_barcode && \\\n";
	print BITMAP "$samtools index -b -@ 24 $outdir/$samp/02.BitBS/$lane_barcode.sort.bam && \\\n";
	print BITMAP "$samtools quickcheck $outdir/$samp/02.BitBS/$lane_barcode.sort.bam && \\\n";
	print BITMAP "echo ==========end at : `date` ========== && \\\n";
	print BITMAP "echo Still_waters_run_deep 1>&2\n";
	close BITMAP;

	open RMCLEAN,">$outdir/$samp/shell/rmCleanReads.$lane_barcode.sh" || die $!;
	print RMCLEAN "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` on `uname -a` ========== && \\\n";
	print RMCLEAN "rm -f $fq1 $fq2 && \\\n";
	print RMCLEAN "echo ==========end at : `date` ========== && \\\n";
	print RMCLEAN "echo Still_waters_run_deep 1>&2\n";
	close RMCLEAN;
	
	`mkdir -p $outdir/$samp/04.QC/BS`;
	open BSMAP,">$outdir/$samp/shell/BS.$lane_barcode.sh" || die $!;
	print BSMAP "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` on `uname -a` ========== && \\\n";
	print BSMAP "export LD_LIBRARY_PATH=/opt/bioinfo/BitMapperBS/v1.0.2.3-centOS7/SSE4.2/htslib_aim/lib && \\\n";
	print BSMAP "$BitMapperBS --search $ref_lambda --seq1 $fq1 --seq2 $fq2 -t $BitCPU --bam --mapstats $outdir/$samp/04.QC/BS/$lane_barcode.BitMapperBS.PE.stats | $samtools addreplacerg -@ $BitCPU -r '\@RG\\tID:$rgid\\tLB:$lib\\tPU:$lane\\tSM:$samp\\tCN:BGI' - | $samtools sort -@ $BitCPU --output-fmt BAM - -o $outdir/$samp/04.QC/BS/$lane_barcode.sort.bam && \\\n";
	print BSMAP "$samtools index -b -@ 8 $outdir/$samp/04.QC/BS/$lane_barcode.sort.bam && \\\n";
	print BSMAP "$samtools quickcheck $outdir/$samp/04.QC/BS/$lane_barcode.sort.bam && \\\n";
	print BSMAP "echo ==========end at : `date` ========== && \\\n";
	print BSMAP "echo Still_waters_run_deep 1>&2\n";
	close BSMAP;
}
close LL;

foreach my $samp (sort keys %samples){
	### Rmdup by sambamba
	my ($instring1, $rmstring1) = ("", "");
	foreach my $lib (sort keys %{$samples{$samp}}){
		foreach my $lane_barcode (sort keys %{$samples{$samp}{$lib}}){
			$instring1 .= " $outdir/$samp/02.BitBS/$lane_barcode.sort.bam ";
			$rmstring1 .= " $outdir/$samp/02.BitBS/$lane_barcode.sort.bam $outdir/$samp/02.BitBS/$lane_barcode.sort.bam.bai $outdir/$samp/02.BitBS/$lane_barcode.sort.bam.csi ";
		}
	}
	my $prefixST = "$outdir/$samp/02.BitBS/$samp";
	open RMDUP,">$outdir/$samp/shell/rmdup.$samp.sh" || die $!;
	print RMDUP "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` on `uname -a` ========== && \\\n";
	print RMDUP "mkdir -p /data/scratch/$samp && \\\n";
    print RMDUP "$sambamba markdup --remove-duplicates -t $sambambaCPU -l 6 --hash-table-size=5000000 --overflow-list-size=20000000 --sort-buffer-size=4096 --tmpdir=/data/scratch/$samp.rmdup.tmp $instring1 $prefixST.rmdup.bam > $prefixST.duplicate_metrics 2>&1 && \\\n";
	print RMDUP "/usr/bin/md5sum $prefixST.rmdup.bam > $prefixST.rmdup.bam.md5sum && \\\n";
	#print RMDUP "rm -rf $rmstring1 && \\\n";
	print RMDUP "rm -rf /data/scratch/$samp.rmdup.tmp && \\\n";
	print RMDUP "echo ==========end at : `date` ========== && \\\n";
	print RMDUP "echo Still_waters_run_deep 1>&2\n";
	close RMDUP;

	### bam qc stat
	open BAMQC,">$outdir/$samp/shell/bamqc.$samp.sh" || die $!;
	print BAMQC "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` on `uname -a` ========== && \\\n";
	print BAMQC "export PATH=$Rbin:\$PATH && \\\n";
	print BAMQC "java -Xmx30g -jar $picard CollectMultipleMetrics R=$reference I=$prefixST.rmdup.bam O=$prefixST PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle PROGRAM=CollectBaseDistributionByCycle PROGRAM=CollectGcBiasMetrics && \\\n";
	print BAMQC "java -Xmx30g -jar $picard CollectWgsMetrics R=$reference I=$prefixST.rmdup.bam O=$prefixST.wgs_metrics READ_LENGTH=$length && \\\n";
	print BAMQC "java -Xmx30g -jar $picard CollectWgsMetrics R=$reference I=$prefixST.rmdup.bam O=$prefixST.CpGisland_metrics READ_LENGTH=$length INTERVALS=$cgibed && \\\n";
	print BAMQC "echo ==========end at : `date` ========== && \\\n";
	print BAMQC "echo Still_waters_run_deep 1>&2\n";
	close BAMQC;

	### Extract CPG, CHG and CHH by MethylDackel
	`mkdir -p $outdir/$samp/03.Methyl`;
	open METHYL,">$outdir/$samp/shell/Methyl.$samp.sh" || die $!;
	print METHYL "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` on `uname -a` ========== && \\\n";
	print METHYL "export LD_LIBRARY_PATH=/opt/bioinfo/samtools/v1.9-centOS6/lib:\$LD_LIBRARY_PATH && \\\n";
	print METHYL "$MethylDackel extract -@ $MethylCPU --CHG --CHH -o $outdir/$samp/03.Methyl/$samp --cytosine_report $reference $outdir/$samp/02.BitBS/$samp.rmdup.bam && \\\n";
	print METHYL "/usr/bin/gzip -f $outdir/$samp/03.Methyl/$samp.cytosine_report.txt && \\\n";
	print METHYL "/usr/bin/gzip -t $outdir/$samp/03.Methyl/$samp.cytosine_report.txt.gz && \\\n";
	print METHYL "/usr/bin/perl $BIN/MethylID.pl $MethylID $outdir/$samp/03.Methyl/$samp.cytosine_report.txt.gz && \\\n";
	print METHYL "/usr/bin/md5sum $outdir/$samp/03.Methyl/$samp.cytosine_report.txt.gz > $outdir/$samp/03.Methyl/$samp.cytosine_report.txt.gz.md5sum && \\\n";
	print METHYL "echo ==========end at : `date` ========== && \\\n";
	print METHYL "echo Still_waters_run_deep 1>&2\n";
	close METHYL;

	### Mbias by MethylDackel
	`mkdir -p $outdir/$samp/04.QC/mbias`;
	open MBIAS,">$outdir/$samp/shell/Mbias.$samp.sh" || die $!;
	print MBIAS "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` on `uname -a` ========== && \\\n";
	print MBIAS "export LD_LIBRARY_PATH=/opt/bioinfo/samtools/v1.9-centOS6/lib:\$LD_LIBRARY_PATH && \\\n";
	print MBIAS "$MethylDackel mbias -@ $MethylCPU --noSVG $reference $outdir/$samp/02.BitBS/$samp.rmdup.bam > $outdir/$samp/04.QC/mbias/$samp.CpG.txt && \\\n";
	print MBIAS "$MethylDackel mbias -@ $MethylCPU --noSVG --noCpG --CHG $reference $outdir/$samp/02.BitBS/$samp.rmdup.bam > $outdir/$samp/04.QC/mbias/$samp.CHG.txt && \\\n";
	print MBIAS "$MethylDackel mbias -@ $MethylCPU --noSVG --noCpG --CHH $reference $outdir/$samp/02.BitBS/$samp.rmdup.bam > $outdir/$samp/04.QC/mbias/$samp.CHH.txt && \\\n";
	print MBIAS "/usr/bin/perl $BIN/Rmbias.pl $outdir/$samp/04.QC/mbias/$samp.CpG.txt $outdir/$samp/04.QC/mbias/$samp.CHG.txt $outdir/$samp/04.QC/mbias/$samp.CHH.txt && \\\n";
	print MBIAS "#Rscript $outdir/$samp/04.QC/mbias/$samp.R && \\\n";
	print MBIAS "rm -f $outdir/$samp/shell/core.* $outdir/$samp/04.QC/mbias/*.txt $outdir/$samp/04.QC/mbias/*.csv $outdir/$samp/04.QC/mbias/*.R && \\\n";
	print MBIAS "echo ==========end at : `date` ========== && \\\n";
	print MBIAS "echo Still_waters_run_deep 1>&2\n";
	close MBIAS;

	open RMBIAS,">$outdir/$samp/04.QC/mbias/$samp.R" || die $!;
	print RMBIAS "#!/path/to/Rscript\n";
	print RMBIAS "CpG_1 <- read.csv(\"$outdir/$samp/04.QC/mbias/$samp.CpG.txt_1.csv\",sep = \"\\t\",header=T)\n";
	print RMBIAS "CpG_2 <- read.csv(\"$outdir/$samp/04.QC/mbias/$samp.CpG.txt_2.csv\",sep = \"\\t\",header=T)\n";
	print RMBIAS "CHG_1 <- read.csv(\"$outdir/$samp/04.QC/mbias/$samp.CHG.txt_1.csv\",sep = \"\\t\",header=T)\n";
	print RMBIAS "CHG_2 <- read.csv(\"$outdir/$samp/04.QC/mbias/$samp.CHG.txt_2.csv\",sep = \"\\t\",header=T)\n";
	print RMBIAS "CHH_1 <- read.csv(\"$outdir/$samp/04.QC/mbias/$samp.CHH.txt_1.csv\",sep = \"\\t\",header=T)\n";
	print RMBIAS "CHH_2 <- read.csv(\"$outdir/$samp/04.QC/mbias/$samp.CHH.txt_2.csv\",sep = \"\\t\",header=T)\n";
	print RMBIAS "out1 <- matrix(ncol=4,nrow=($length-$mbias1))\n";
	print RMBIAS "out2 <- matrix(ncol=4,nrow=($length-$mbias2))\n";
	print RMBIAS "for (i in 1:($length-$mbias1)) {\n";
	print RMBIAS "	position=i+$mbias1;\n";
	print RMBIAS "	CpG_methylation_read1=(CpG_1[i,4]+CpG_1[(2*$length+1-2*$mbias1-i),4])/0.02\n";
	print RMBIAS "	CHG_methylation_read1=(CHG_1[i,4]+CHG_1[(2*$length+1-2*$mbias1-i),4])/0.02\n";
	print RMBIAS "	CHH_methylation_read1=(CHH_1[i,4]+CHH_1[(2*$length+1-2*$mbias1-i),4])/0.02\n";
	print RMBIAS "	out1[i,] <- cbind(position,CpG_methylation_read1,CHG_methylation_read1,CHH_methylation_read1)\n}\n";
	print RMBIAS "for (i in 1:($length-$mbias2)) {\n";
	print RMBIAS "	position=i+$mbias2;\n";
	print RMBIAS "	CpG_methylation_read2=(CpG_2[($length-$mbias2+i),4]+CpG_2[($length+1-$mbias2-i),4])/0.02\n";
	print RMBIAS "	CHG_methylation_read2=(CHG_2[($length-$mbias2+i),4]+CHG_2[($length+1-$mbias2-i),4])/0.02\n";
	print RMBIAS "	CHH_methylation_read2=(CHH_2[($length-$mbias2+i),4]+CHH_2[($length+1-$mbias2-i),4])/0.02\n";
	print RMBIAS "	out2[i,] <- cbind(position,CpG_methylation_read2,CHG_methylation_read2,CHH_methylation_read2)\n";
	print RMBIAS "}\ncolnames(out1)<-c(\"position\",\"CpG_methylation_read1\",\"CHG_methylation_read1\",\"CHH_methylation_read1\")\n";
	print RMBIAS "colnames(out2)<-c(\"position\",\"CpG_methylation_read2\",\"CHG_methylation_read2\",\"CHH_methylation_read2\")\n";
	print RMBIAS "write.csv(out1,file=\"$outdir/$samp/04.QC/mbias/$samp\_mbias1.csv\",quote = F,row.names = F)\n";
	print RMBIAS "write.csv(out2,file=\"$outdir/$samp/04.QC/mbias/$samp\_mbias2.csv\",quote = F,row.names = F)\n";
	print RMBIAS "mbias1 <- read.csv(\"$outdir/$samp/04.QC/mbias/$samp\_mbias1.csv\",header = T)\n";
	print RMBIAS "mbias2 <- read.csv(\"$outdir/$samp/04.QC/mbias/$samp\_mbias2.csv\",header = T)\n";
	print RMBIAS "p1=spline(mbias1\$position,mbias1\$CpG_methylation_read1,n=1000)\n";
	print RMBIAS "p2=spline(mbias1\$position,mbias1\$CHG_methylation_read1,n=1000)\n";
	print RMBIAS "p3=spline(mbias1\$position,mbias1\$CHH_methylation_read1,n=1000)\n";
	print RMBIAS "p4=spline(mbias2\$position,mbias2\$CpG_methylation_read2,n=1000)\n";
	print RMBIAS "p5=spline(mbias2\$position,mbias2\$CHG_methylation_read2,n=1000)\n";
	print RMBIAS "p6=spline(mbias2\$position,mbias2\$CHH_methylation_read2,n=1000)\n";
	print RMBIAS "png(filename=\"$outdir/$samp/04.QC/mbias/$samp\_mbias.png\",width = 850,height = 428)\n";
	print RMBIAS "par(mai=c(1,1,0.5,3))\n";
	print RMBIAS "plot(p1,col=\"royalblue4\",type=\"l\",xlim=c(1,$length),ylim=c(0,100),lwd=2,xlab=\"position in reads (bp)\",ylab=\"methylation rate (%)\",main=\"$samp M-bias\",col.main=\"Black\",font.main=2,lty=1)\n";
	print RMBIAS "lines(p2,col=\"orange2\",type=\"l\",xlim=c(1,$length),ylim=c(0,1),lwd=2,xlab=\"position in reads (bp)\",ylab=\"methylation rate (%)\",main=\"$samp M-bias\",col.main=\"Black\",font.main=2,lty=1)\n";
	print RMBIAS "lines(p3,col=\"cyan\",type=\"l\",xlim=c(1,$length),ylim=c(0,1),lwd=2,xlab=\"position in reads (bp)\",ylab=\"methylation rate (%)\",main=\"$samp M-bias\",col.main=\"Black\",font.main=2,lty=1)\n";
	print RMBIAS "lines(p4,col=\"royalblue4\",type=\"l\",xlim=c(1,$length),ylim=c(0,1),lwd=2,xlab=\"position in reads (bp)\",ylab=\"methylation rate (%)\",main=\"13-219A M-bias\",col.main=\"Black\",font.main=2,lty=2)\n";
	print RMBIAS "lines(p5,col=\"orange2\",type=\"l\",xlim=c(1,$length),ylim=c(0,1),lwd=2,xlab=\"position in reads (bp)\",ylab=\"methylation rate (%)\",main=\"$samp M-bias\",col.main=\"Black\",font.main=2,lty=2)\n";
	print RMBIAS "lines(p6,col=\"cyan\",type=\"l\",xlim=c(1,$length),ylim=c(0,1),lwd=2,xlab=\"position in reads (bp)\",ylab=\"methylation rate (%)\",main=\"$samp M-bias\",col.main=\"Black\",font.main=2,lty=2)\n";
	print RMBIAS "xy <- par(\"usr\")\n";
	print RMBIAS "legend(x=xy[2]-xinch(-0.1),y=xy[3]-yinch(-2.5),legend=c(\"CpG_methylation_read1\",\"CHG_methylation_read1\",\"CHH_methylation_read1\",\"CpG_methylation_read2\",\"CHG_methylation_read2\",\"CHH_methylation_read2\"),col=c(\"royalblue4\",\"orange2\",\"cyan\",\"royalblue4\",\"orange2\",\"cyan\"),lwd=0.5,lty=c(1,1,1,2,2,2),bty=\"n\",cex = 1,xpd = TRUE)\n";
	print RMBIAS "dev.off()\n";
	close RMBIAS;

	### Calculated conversion efficiency
	my ($instring2, $rmstring2) = ("", "");
	foreach my $lib (sort keys %{$samples{$samp}}){
		foreach my $lane_barcode (sort keys %{$samples{$samp}{$lib}}){
			$instring2 .= " $outdir/$samp/04.QC/BS/$lane_barcode.sort.bam ";
			$rmstring2 .= " $outdir/$samp/04.QC/BS/$lane_barcode.sort.bam $outdir/$samp/04.QC/BS/$lane_barcode.sort.bam.bai $outdir/$samp/04.QC/BS/$lane_barcode.sort.bam.csi ";
		}
	}
	open BSRMDUP,">$outdir/$samp/shell/BS.rmdup.$samp.sh" || die $!;
	print BSRMDUP "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` on `uname -a` ========== && \\\n";
	print BSRMDUP "export LD_LIBRARY_PATH=/opt/bioinfo/samtools/v1.9-centOS6/lib:\$LD_LIBRARY_PATH && \\\n";
	print BSRMDUP "$sambamba markdup --remove-duplicates -t $sambambaCPU -l 6 --hash-table-size=1000000 --overflow-list-size=2000000 --sort-buffer-size=4096 --tmpdir=$outdir/$samp/04.QC/BS/BS.$samp.rmdup.tmp $instring2 $outdir/$samp/04.QC/BS/BS.$samp.rmdup.bam > $outdir/$samp/04.QC/BS/BS.$samp.duplicate_metrics 2>&1 && \\\n";
	print BSRMDUP "$MethylDackel extract --CHG --CHH -o $outdir/$samp/04.QC/BS/BS.$samp --cytosine_report $ref_lambda $outdir/$samp/04.QC/BS/BS.$samp.rmdup.bam && \\\n";
	print BSRMDUP "awk 'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{mc+=\$4;umc+=\$5}END{print \"conversion efficiency\",umc/(mc+umc)}' $outdir/$samp/04.QC/BS/BS.$samp.cytosine_report.txt > $outdir/$samp/02.BitBS/conversion.stat && \\\n";
	print BSRMDUP "rm -rf $rmstring2 && \\\n";
	print BSRMDUP "rm -rf $outdir/$samp/04.QC/BS/BS.$samp.rmdup.tmp && \\\n";
	print BSRMDUP "echo ==========end at : `date` ========== && \\\n";
	print BSRMDUP "echo Still_waters_run_deep 1>&2\n";
	close BSRMDUP;

	### reads stat
	`mkdir -p $outdir/shell_run`;
	open ALLSTATSH, ">$outdir/shell_run/stat.$subProjectName.sh" || die $!;
	print ALLSTATSH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` on `uname -a` ========== && \\\n";
	print ALLSTATSH "/usr/bin/perl $Barcode_stat -list $List -outdir $outdir > $outdir/$subProjectName.Barcode_stat.txt && \\\n";
	print ALLSTATSH "/usr/bin/perl $Sample_stat -list $List -outdir $outdir -barcodestat $outdir/$subProjectName.Barcode_stat.txt > $outdir/$subProjectName.Sample_stat.txt && \\\n";
	print ALLSTATSH "echo ==========end at : `date` ========== && \\\n";
	print ALLSTATSH "echo Still_waters_run_deep 1>&2\n";
	close ALLSTATSH;
}
	
##################################################moniotor########################################################
open EDGE,">$outdir/shell_run/edge.$subProjectName.list" or die$!;

for my $samp (sort keys %samples){
	if ($Use_SOAPnuke=~/true/i){
		for my $lib(sort keys %{$samples{$samp}}){
			for my $lane_barcode (sort keys %{$samples{$samp}{$lib}}){
				print EDGE "$outdir/$samp/shell/Clean.$lane_barcode.sh:6G:6cpu\t$outdir/$samp/shell/fqcheck.$lane_barcode.sh:3G:1cpu\n";
				print EDGE "$outdir/$samp/shell/Clean.$lane_barcode.sh:6G:6cpu\t$outdir/$samp/shell/Alignment.$lane_barcode.sh:20G:$BitCPU"."cpu\n";
				print EDGE "$outdir/$samp/shell/Clean.$lane_barcode.sh:6G:6cpu\t$outdir/$samp/shell/BS.$lane_barcode.sh:20G:$BitCPU"."cpu\n";
				if ($rmCleanReads=~/true/i){
					print EDGE "$outdir/$samp/shell/fqcheck.$lane_barcode.sh:3G:1cpu\t$outdir/$samp/shell/rmCleanReads.$lane_barcode.sh:1G:1cpu\n";
					print EDGE "$outdir/$samp/shell/rmdup.$samp.sh:30G:$sambambaCPU"."cpu\t$outdir/$samp/shell/rmCleanReads.$lane_barcode.sh:1G:1cpu\n";
					print EDGE "$outdir/$samp/shell/BS.rmdup.$samp.sh:3G:1cpu\t$outdir/$samp/shell/rmCleanReads.$lane_barcode.sh:1G:1cpu\n";
				}
			}
		}
	}	
	for my $lib(sort keys %{$samples{$samp}}){
		for my $lane_barcode (sort keys %{$samples{$samp}{$lib}}){
			print EDGE "$outdir/$samp/shell/Alignment.$lane_barcode.sh:20G:$BitCPU"."cpu\t$outdir/$samp/shell/rmdup.$samp.sh:30G:$sambambaCPU"."cpu\n";
			print EDGE "$outdir/$samp/shell/BS.$lane_barcode.sh:20G:$BitCPU"."cpu\t$outdir/$samp/shell/BS.rmdup.$samp.sh:3G:1cpu\n";
		}
	}
	print EDGE "$outdir/$samp/shell/rmdup.$samp.sh:30G:$sambambaCPU"."cpu\t$outdir/$samp/shell/Methyl.$samp.sh:3G:$MethylCPU"."cpu\n";
	print EDGE "$outdir/$samp/shell/rmdup.$samp.sh:30G:$sambambaCPU"."cpu\t$outdir/$samp/shell/Mbias.$samp.sh:3G:$MethylCPU"."cpu\n";
	print EDGE "$outdir/$samp/shell/rmdup.$samp.sh:30G:$sambambaCPU"."cpu\t$outdir/$samp/shell/bamqc.$samp.sh:30G:1cpu\n";
	print EDGE "$outdir/$samp/shell/bamqc.$samp.sh:30G:1cpu\t$outdir/shell_run/stat.$subProjectName.sh:1G:1cpu\n";
}
close EDGE;

open RUNSH,">$outdir/shell_run/run.$subProjectName.sh" or die$!;
print RUNSH "$monitor taskmonitor -q $queue -P $priority -p $subProjectName -f 0 -s Still_waters_run_deep -i $outdir/shell_run/edge.$subProjectName.list\n";
close RUNSH;
