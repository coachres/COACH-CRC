#!perl -w
use strict;
use Getopt::Long;
my $usage=<<'USAGE';

Somatic mutation calling pipeline of hg38 for Sentieon Platform.
Version 1.0: 2022-03-20
Author : luoshuzhen, CANCER INSTITUTE

Options:
    -list            STR  bam file list.  <require>
    -pair            STR  pair sample info.  <require>
    -outdir          STR  the output directory.  <require>
    -rundir          STR  the monitor directory.  <require>
    -proj            STR  pair sample info.  <require>
    -local           STR  copy data to local directory.  [default:F]
    -h                    Help Information

USAGE

my ($fileList,$pairList,$somaticdir,$runshelldir,$projectName,$local)=@ARGV;

GetOptions(
    'list=s' => \$fileList,
    'pair=s' => \$pairList,
    'outdir=s' => \$somaticdir,
    'rundir=s' => \$runshelldir,
    'proj=s' => \$projectName,
    'local=s' => \$local,
    'h|?' => sub{die "$usage\n";},
);

die "$usage\n" unless (defined($fileList) && defined($pairList) && defined($somaticdir) && defined($runshelldir) && defined($projectName));

$local ||= "F";

#========================== tools ==========================
my $Bin_CentOS7         = "/opt/bioinfo/tools";
my $bgzip               = "$Bin_CentOS7/samtools/v1.12/bin/bgzip";
my $tabix               = "$Bin_CentOS7/samtools/v1.12/bin/tabix";
my $bcftools            = "$Bin_CentOS7/samtools/v1.12/bin/bcftools";
my $bedtools            = "/opt/bioinfo/tools/bedtools/v2.27.1/bin/bedtools";
my $Vt                  = "/opt/bioinfo/tools/vt/vt-0.57721";
my $expr_conf           = "$Bin_CentOS7/multTools_somatic_snvs_indels_ensemble.bcftools_expr.conf";
my $msisensor_pro       = "/opt/bioinfo/tools/MSIsensor-pro/v1.2.0"; #
my $MuSE                = "/opt/bioinfo/tools/MuSE/v1.0/MuSE";
my $vcfutils            = "$Bin_CentOS7/samtools/v1.12/bin/vcfutils.pl";
my $java8               = "/opt/bioinfo/tools/java/jdk1.8.0_202/bin/java";
my $GATK4               = "/opt/bioinfo/tools/gatk-4.1.0.0/gatk";
my $GATK3               = "/opt/bioinfo/tools/gatk-4.1.0.0/gatk";
my $bam_readcount       = "/opt/bioinfo/tools/bam-readcount/v0.8.0/bin/bam-readcount";
my $bam_somaticsniper   = "/opt/bioinfo/tools/somatic-sniper/v1.0.5";
my $Manta               = "/opt/bioinfo/tools/Manta/v1.5.0/bin/configManta.py";
my $Strelka2            = "/opt/bioinfo/tools/strelka/v2.9.10/bin/configureStrelkaSomaticWorkflow.py";
my $Lancet              = "/opt/bioinfo/tools/lancet/v1.1.0/lancet";
my $Svaba               = "/opt/bioinfo/tools/SvABA/v0.2.1/bin/svaba";
my $bammathPATH         = "/opt/bioinfo/tools/tools/BAM-matcher/v20170604";
#my $pythonBin           = "/opt/python/bin";

##========================== database ==========================##
my $HG38_PATH           = "/opt/bioinfo/database/gatk_bundle/hg38";
my $hapmapAF            = "$HG38_PATH/hapmap_3.3_grch38_pop_stratified_af.vcf.gz";
my $reffa               = "$HG38_PATH/Homo_sapiens_assembly38.fasta";
my $refdict             = "$HG38_PATH/Homo_sapiens_assembly38.dict";
my $reffai              = "$HG38_PATH/Homo_sapiens_assembly38.fasta.fai";
my $cosmic              = "$HG38_PATH/cosmic_v92.hg38.vcf.gz";
my $dbsnp               = "$HG38_PATH/dbsnp_151.hg38.vcf.gz";
my $dbsnpcommon         = "$HG38_PATH/dbsnp_151.common.hg38.vcf.gz";
my $gnomadGATK          = "$HG38_PATH/af-only-gnomad.hg38.vcf.gz";
my $svdprefix           = "/VerifyBamID/resource/1000g.phase3.100k.b38.vcf.gz.dat";   # SVD related files prefix for VerifyBamID2.
my $pcgr_data           = "/opt/bioinfo/tools/PCGR/pcgr-0.9.1/data"; # PCGR v0.9.1, data-20201123
my $CGC_bed             = "/opt/bioinfo/database/hg38/Cosmic/v94/cancer_gene_census.bed"; # cancer_gene_census, Cosmic v94
my $CancerIT_hts_ref    = "/opt/bioinfo/tools/CancerIT/hts-ref"; # for cram reading, REF_PATH and REF_CACHE

##========================== genome region for mutations ====================##
my $genomebed = "/opt/bioinfo/04.Somatic/Region/chrmosomes.bed";
my $gatkintervals = "/opt/bioinfo/04.Somatic/Region/chrmosomes.bed.intervals";
my $msisite = "/opt/bioinfo/tools/MSIsensor-pro/v1.2.0/Homo_sapiens_assembly38.microsate.sites";

##========================== input and out path ==========================##

unless(-e $runshelldir){system("mkdir -p $runshelldir")}

##========================== MONITOR ==========================##

my %filelist;
open IA, $fileList || die $!;
while(my $L1=<IA>){
    chomp($L1);
    my ($name,$file) = split/\t/,$L1;
    $filelist{$name} = $file;
}
close IA;

open PA, $pairList || die $!;
while(my $L2=<PA>){
    chomp($L2);
    my ($nname,$tname) = split/\t/,$L2;
    if(exists $filelist{$nname} && exists $filelist{$tname}){
        my $pairname = "$nname-VS-$tname";
        my $pairDir = "$somaticdir/$pairname";
        my $pairshdir = "$somaticdir/$pairname/shell";
        unless(-e $pairDir){system("mkdir -p $pairDir")}
        unless(-e $pairshdir){system("mkdir -p $pairshdir")}
        my $normalbam = $filelist{$nname};
        my $tumorbam = $filelist{$tname};
        my $nfilename = (split/\//,$normalbam)[-1];
        my $tfilename = (split/\//,$tumorbam)[-1];
        my ($normalbamLocal,$tumorbamLocal,$pairDirLocal,$combinaSH);
        my $LocalDir = "/data/scratch/";
        if($local eq "T"){
            $pairDirLocal = "$LocalDir/$pairname";
            $normalbamLocal = "$pairDirLocal/bamFile/$nfilename";
            $tumorbamLocal = "$pairDirLocal/bamFile/$tfilename";
             }else{
                $normalbamLocal = $normalbam;
                $tumorbamLocal = $tumorbam;
                $pairDirLocal = $pairDir;
            }
        ## somatic mutation calling
        &MuSECalling($pairDir,$pairname,$normalbam,$tumorbam,$pairshdir,$genomebed);
        &MuTect2Calling($pairDir,$pairname,$nname,$tname,$normalbam,$tumorbam,$pairshdir,$gatkintervals);
        &SomaticSniperCalling($pairDir,$pairname,$nname,$tname,$normalbam,$tumorbam,$pairshdir,$genomebed);

        my $step1 = "echo ==========Manta Calling start at : `date` ========== && \\\n";
    $step1 .= &MantaCalling($pairDirLocal,$pairDir,$pairname,$normalbamLocal,$tumorbamLocal,$pairshdir,$genomebed.".gz");
    if( -e "$pairDir/Manta/results/variants/candidateSmallIndels.vcf.gz"){
    }else{
         $step1 .= "rsync -atvP $pairDirLocal/Manta $pairDir && \\\n";
    }

    my $step2 = "echo ==========Strelka2 Calling start at : `date` ==========&& \\\n";
    $step2 .= &Strelka2Calling($pairDirLocal,$pairname,$normalbamLocal,$tumorbamLocal,$pairshdir,$genomebed.".gz","true");
    $step2 .= "rsync -atvP $pairDirLocal/Strelka2 $pairDir && \\\n";

    my $step3 = "echo ==========Svaba Calling start at : `date` ==========&& \\\n";
    $step3 .= &runSvaba($pairDirLocal,$pairname,$normalbamLocal,$tumorbamLocal,$pairshdir,$genomebed);
    $step3 .= "rsync -atvP $pairDirLocal/Svaba $pairDir && \\\n";

    my $step4 = "echo ==========Lancet Calling start at : `date` ==========&& \\\n";
    $step4 .= &LancetCalling($pairDirLocal,$pairDir,$pairname,$normalbamLocal,$tumorbamLocal,$pairshdir,$genomebed);
    $step4 .= "rsync -atvP $pairDirLocal/Lancet $pairDir && \\\n";

    my $step5 = "echo ==========MSI Calling start at : `date` ==========\n";
    $step5 .= &MSIsensorPro($pairDirLocal,$normalbamLocal,$tumorbamLocal,$pairshdir,$tname,$msisite);
    $step5 .= "rsync -atvP $pairDirLocal/MSIsensorPro $pairDir && \\\n";

    open COM, ">$pairshdir/combinaTools.$pairname.sh" || die $!;

    print COM "#!/bin/bash\n";
    print COM "set -euo pipefail\n";
    if($local eq "T"){
        print COM "echo ==========Copy bam files start at : `date` on `hostname` ==========\n";
        print COM "mkdir -p $pairDirLocal/bamFile && \\\n";
        print COM "rm -rf $pairDirLocal/bamFile/* && \\\n";
        print COM "cp -i $normalbam $pairDirLocal/bamFile/ && \\\n";
        print COM "cp -i $normalbam.bai $pairDirLocal/bamFile/ && \\\n";
        print COM "cp -i $tumorbam $pairDirLocal/bamFile/ && \\\n";
        print COM "cp -i $tumorbam.bai $pairDirLocal/bamFile/ && \\\n";
        print COM "echo ==========Copy bam files end at : `date` ========== && \\\n";
    }else{ }
    print COM $step1;
    print COM $step2;
    print COM $step3;
    print COM $step4;
    #print COM $step5;
    print COM "rm -rf $pairDirLocal/bamFile/*.bam && echo -e bam files had removed\n";
    print COM "if [ ! -e \"$pairDirLocal/Lancet/$pairname.Lancet.snv_shortindel.vcf.gz\" ];then  echo Error at `date` on `hostname` >$pairshdir/eror.log && rm -rf $pairDirLocal/bamFile/*.bam;fi\n";
    print COM "echo ========== all were end at : `date` on `hostname` ========== \n";

    ## sample QC
    #&somaticQC($pairDir,$pairname,$normalbam,$tumorbam,$pairshdir);

    ## MSI detect
    #&MSIsensorPro($pairDir,$normalbam,$tumorbam,$pairshdir,$tname,$msisite);

}
close PA;

#========================== MSI detection by MSIsensor-pro ==========================##
## base on https://github.com/xjtu-omics/msisensor-pro
sub MSIsensorPro{
    my ($pairdir, $nbam, $tbam,$shdir, $name, $microsatellites) = @_;
    my $dir = "$pairdir/MSIsensorPro";
    my $shell = "$shdir/MSIsensorPro.$name.sh";
    #unless(-e $dir){system("mkdir -p $dir")}
    #open SH, ">$shell" || die $!;
    #print SH &echoSTARTstring($shell);
    my $msiSH = "mkdir -p $dir && \\\n";
    $msiSH .= "$msisensor_pro/msisensor-pro msi -d $microsatellites -n $nbam -t $tbam -b 16 -c 15 -o $dir/$name.result && \\\n";
    $msiSH .= "gzip -f $dir/$name.result_dis $dir/$name.result_somatic && \\\n";
    #$msiSH .= &echoENDstring($shell);
    if($local eq "T"){
    }else{
        open SH, ">$shell" || die $!;
        print SH &echoSTARTstring($shell);
        print SH $msiSH;
        print SH &echoENDstring($shell);
        close SH;
    }
    return($msiSH);
}

##========================== calling somatic snv by MuSE ============================##
## tumor.bam matched_normal.bam
## filter raw snv.  From experience, we suggest using calls up to Tier 5 for WGS data.
sub MuSECalling{
    my ($pairdir,$samplePair,$nbam,$tbam,$shdir,$bed)=@_;
    my $shell01 = "$shdir/MuSE.$samplePair.sh";
    my $dir = "$pairdir/MuSE";
    unless(-e $dir){system("mkdir -p $dir")}
    open SH1, ">$shell01" || die $!;
    print SH1 &echoSTARTstring($shell01);
    print SH1 "$MuSE call -f $reffa -l $bed -O $dir/$samplePair.raw $tbam $nbam && \\\n";
    print SH1 "colnum=`cat $dir/$samplePair.raw.MuSE.txt |tail -1|awk -F \"\\t\" '{print NF}'`\n";
    print SH1 "$MuSE sump -I $dir/$samplePair.raw.MuSE.txt -O $dir/$samplePair.raw.MuSE.vcf -D $dbsnp -G && \\\n";
    print SH1 "$bgzip -f $dir/$samplePair.raw.MuSE.txt && \\\n";
    print SH1 "$bgzip -f $dir/$samplePair.raw.MuSE.vcf && \\\n";
    print SH1 "$tabix -f -p vcf $dir/$samplePair.raw.MuSE.vcf.gz && \\\n";
    print SH1 "zcat $dir/$samplePair.raw.MuSE.vcf.gz | awk '\$1~/^#/ || \$7~/PASS|Tier1|Tier2|Tier3|Tier4|Tier5/ {print}' | $bgzip > $dir/$samplePair.MuSE.snv.vcf.gz && \\\n";
    print SH1 "$tabix -f -p vcf $dir/$samplePair.MuSE.snv.vcf.gz && \\\n";
    print SH1 &getBed("$dir/$samplePair.MuSE.snv.vcf.gz","$dir/$samplePair.MuSE.called.bed");
    print SH1 &echoENDstring($shell01);
    close SH1;
}

##========================== calling somatic snv by MuTect2 ============================##
sub MuTect2Calling{
    my ($pairdir,$samplePair,$normal,$tumor,$nbam,$tbam,$shdir,$bed)=@_;
    my $shell02 = "$shdir/MuTect2.$samplePair.sh";
    my $MuTect2Dir = "$pairdir/MuTect2";
    unless(-e $MuTect2Dir){system("mkdir -p $MuTect2Dir");system("mkdir -p $MuTect2Dir/tmp")}
    open SH2, ">$shell02" || die $!;
    print SH2 &echoSTARTstring($shell02);
    print SH2 "export PATH=/opt/bioinfo/tools/java/jdk1.8.0_202/bin:\$PATH && \\\n";
    print SH2 "$GATK4 --java-options '-Xmx6g -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$MuTect2Dir/tmp' Mutect2 -R $reffa -I $tbam -tumor $tumor -I $nbam -normal $normal --germline-resource $gnomadGATK --af-of-alleles-not-in-resource 3.125e-05 -L $bed -O $MuTect2Dir/$samplePair.MuTect2.raw.vcf.gz && \\\n";
    print SH2 "$GATK4 --java-options '-Xmx6g -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$MuTect2Dir/tmp' FilterMutectCalls -V $MuTect2Dir/$samplePair.MuTect2.raw.vcf.gz -O $MuTect2Dir/$samplePair.MuTect2.Filter.vcf.gz && \\\n";
    print SH2 "zcat $MuTect2Dir/$samplePair.MuTect2.Filter.vcf.gz|awk '\$1~/^#/ || \$7~/PASS/ {print}' | $bgzip > $MuTect2Dir/$samplePair.MuTect2.snv_shortindel.vcf.gz && \\\n";
    print SH2 "$tabix -f -p vcf $MuTect2Dir/$samplePair.MuTect2.snv_shortindel.vcf.gz && \\\n";
    print SH2 "rm -rf $MuTect2Dir/tmp && \\\n";
    print SH2 &getBed("$MuTect2Dir/$samplePair.MuTect2.snv_shortindel.vcf.gz","$MuTect2Dir/$samplePair.MuTect2.called.bed");
    print SH2 &echoENDstring($shell02);
    close SH2;
}

##========================== calling somatic snv by SomaticSniper ============================##
sub SomaticSniperCalling{
    my ($pairdir,$samplePair,$normal,$tumor,$nbam,$tbam,$shdir,$bed)=@_;
    my $SomaticSniperDir = "$pairdir/SomaticSniper";
    unless(-e $SomaticSniperDir){system("mkdir -p $SomaticSniperDir")}
    my $shell03 = "$shdir/SomaticSniper.$samplePair.sh";
    open SH3, ">$shell03" || die $!;
    print SH3 &echoSTARTstring($shell03);
    print SH3 "$bam_somaticsniper/bin/bam-somaticsniper -q 1 -Q 40 -L -G -F vcf -t $tumor -n $normal -f $reffa $tbam $nbam $SomaticSniperDir/$samplePair.SomaticSniper.raw.vcf.gz && \\\n";
    print SH3 "$bcftools mpileup -A -B -f $reffa $nbam -R $bed|$bcftools call -c|grep -v \"#\"|$vcfutils varFilter -Q 20 > $SomaticSniperDir/$normal.indel.pileup && \\\n";
    print SH3 "$bcftools mpileup -A -B -f $reffa $tbam -R $bed|$bcftools call -c|grep -v \"#\"|$vcfutils varFilter -Q 20 > $SomaticSniperDir/$tumor.indel.pileup && \\\n";
    print SH3 "perl $bam_somaticsniper/src/scripts/snpfilter.pl --snp-file $SomaticSniperDir/$samplePair.SomaticSniper.raw.vcf.gz --indel-file $SomaticSniperDir/$normal.indel.pileup --out-file $SomaticSniperDir/$normal.SNPfilter && \\\n";
    print SH3 "perl $bam_somaticsniper/src/scripts/snpfilter.pl --snp-file $SomaticSniperDir/$normal.SNPfilter --indel-file $SomaticSniperDir/$tumor.indel.pileup --out-file $SomaticSniperDir/$samplePair.SNPfilter && \\\n";
    print SH3 "perl $bam_somaticsniper/src/scripts/prepare_for_readcount.pl --snp-file $SomaticSniperDir/$samplePair.SNPfilter --out-file $SomaticSniperDir/$samplePair.SNPfilter.pos && \\\n";
    print SH3 "$bam_readcount -w 1 -b 15 -q 1 -f $reffa -l $SomaticSniperDir/$samplePair.SNPfilter.pos $tbam > $SomaticSniperDir/$samplePair.readcounts && \\\n";
    print SH3 "perl $bam_somaticsniper/src/scripts/fpfilter.pl --snp-file $SomaticSniperDir/$samplePair.SNPfilter --readcount-file $SomaticSniperDir/$samplePair.readcounts && \\\n";
    print SH3 "perl $bam_somaticsniper/src/scripts/highconfidence.pl --min-mapping-quality 40 --min-somatic-score 40 --snp-file $SomaticSniperDir/$samplePair.SNPfilter.fp_pass --lq-output $SomaticSniperDir/$samplePair.low_qual.vcf.gz --out-file $SomaticSniperDir/$samplePair.SomaticSniper.snv.vcf && \\\n";
    print SH3 "$bgzip -f $SomaticSniperDir/$samplePair.SomaticSniper.snv.vcf && \\\n";
    print SH3 "$tabix -f -p vcf $SomaticSniperDir/$samplePair.SomaticSniper.snv.vcf.gz && \\\n";
    print SH3 &getBed("$SomaticSniperDir/$samplePair.SomaticSniper.snv.vcf.gz","$SomaticSniperDir/$samplePair.SomaticSniper.called.bed");
    print SH3 &echoENDstring($shell03);
    close SH3;
}

##====================== calling structural variants (SVs) and indel by Manta ===============##
sub MantaCalling{
    my ($pairdir,$pairdirshare,$samplePair,$nbam,$tbam,$shdir,$bed)=@_;
    my $MantaDir= "$pairdir/Manta";
    #unless(-e $MantaDir){system("mkdir -p $MantaDir")}
    my $shell04 = "$shdir/Manta.$samplePair.sh";
    #open SH4, ">$shell04" || die $!;
    #print SH4 &echoSTARTstring($shell04);
    my $MantaSH = "mkdir -p $MantaDir && \\\n";
    $MantaSH .= "rm -rf $MantaDir/* && \\\n";
    $MantaSH .= "/usr/bin/python $Manta --normalBam $nbam --tumorBam $tbam --referenceFasta $reffa --runDir $MantaDir --callRegions $bed --generateEvidenceBam --outputContig && \\\n";
    $MantaSH .= "/usr/bin/python $MantaDir/runWorkflow.py -m local -j 32 && \\\n";
    $MantaSH .= "rm -rf $MantaDir/workspace && \\\n";
    #$MantaSH .= &echoENDstring($shell04);
    if($local eq "T"){
    }else{
        open SH4, ">$shell04" || die $!;
        print SH4 &echoSTARTstring($shell04);
        print SH4 $MantaSH;
        print SH4 &echoENDstring($shell04);
        close SH4;
    }
    if( -e "$pairdirshare/Manta/results/variants/candidateSmallIndels.vcf.gz"){
        $MantaSH = "mkdir -p $MantaDir/results/variants && \\\n";
        $MantaSH .= "cp -r $pairdirshare/Manta/results/variants/candidateSmallIndels.vcf.gz $MantaDir/results/variants/candidateSmallIndels.vcf.gz && \\\n";
        $MantaSH .= "cp -r $pairdirshare/Manta/results/variants/candidateSmallIndels.vcf.gz.tbi $MantaDir/results/variants/candidateSmallIndels.vcf.gz.tbi && \\\n";
    }else{}
    return($MantaSH);
}

##========================== calling somatic snv by Strelka2 ============================##
sub Strelka2Calling{
    my ($pairdir,$samplePair,$nbam,$tbam,$shdir,$bed,$Use_MantaAndStrelka2)=@_;
    my $Strelka2Dir = "$pairdir/Strelka2";
    #unless(-e $Strelka2Dir){system("mkdir -p $Strelka2Dir")}
    my $shell05 = "$shdir/Strelka2.$samplePair.sh";
    #open SH5, ">$shell05" || die $!;
    #print SH5 &echoSTARTstring($shell05);
    my $Strelka2SH = "mkdir -p $Strelka2Dir && \\\n";
    $Strelka2SH .= "rm -fr $Strelka2Dir/* && \\\n";
    my $add_para = "";
    if ($Use_MantaAndStrelka2=~/true/i){
        $add_para = "--indelCandidates $pairdir/Manta/results/variants/candidateSmallIndels.vcf.gz";
    }else{
    }
    $Strelka2SH .= "/usr/bin/python $Strelka2 --normalBam $nbam --tumorBam $tbam --ref $reffa --callRegions $bed $add_para --runDir $Strelka2Dir && \\\n";
    $Strelka2SH .= "/usr/bin/python $Strelka2Dir/runWorkflow.py -m local -j 32 && \\\n";
    $Strelka2SH .= "zcat $Strelka2Dir/results/variants/somatic.snvs.vcf.gz | awk '\$1~/^#/ || \$7~/PASS/ {print}' | $bgzip > $Strelka2Dir/$samplePair.Strelka2.snv.vcf.gz && \\\n";
    $Strelka2SH .= "zcat $Strelka2Dir/results/variants/somatic.indels.vcf.gz | awk '\$1~/^#/ || \$7~/PASS/ {print}' | $bgzip > $Strelka2Dir/$samplePair.Strelka2.indel.vcf.gz && \\\n";
    $Strelka2SH .= &getBed("$Strelka2Dir/$samplePair.Strelka2.snv.vcf.gz","$Strelka2Dir/$samplePair.Strelka2.snv.called.bed");
    $Strelka2SH .= &getBed("$Strelka2Dir/$samplePair.Strelka2.indel.vcf.gz","$Strelka2Dir/$samplePair.Strelka2.indel.called.bed");
    $Strelka2SH .= "rm -rf $Strelka2Dir/workspace && \\\n";
    #$Strelka2SH .= &echoENDstring($shell05);
    if($local eq "T"){
    }else{
        open SH5, ">$shell05" || die $!;
        print SH5 &echoSTARTstring($shell05);
        print SH5 $Strelka2SH;
        print SH5 &echoENDstring($shell05);
        close SH5;
    }
    return($Strelka2SH);
}

##====================== calling somatic snv/indel by Lancet =========================##
sub LancetCalling{
    my ($pairdir,$pairdirShare,$samplePair,$nbam,$tbam,$shdir,$bed)=@_;
    my $LancetDir = "$pairdir/Lancet";
    #unless(-e $LancetDir){system("mkdir -p $LancetDir")}
    my $shell06 = "$shdir/Lancet.$samplePair.sh";
    #open SH6, ">$shell06" || die $!;
    #print SH6 &echoSTARTstring($shell06);
    my $GetBedSH = "mkdir -p $LancetDir && \\\n";
    $GetBedSH .= "rm -fr $LancetDir/$samplePair.NonLancet.called.bed && \\\n";
    foreach my $soft('Strelka2','MuSE','MuTect2','SomaticSniper','Svaba'){
            if($soft eq "Strelka2"){
                $GetBedSH .= "cat $pairdir/$soft/$samplePair.Strelka2.snv.called.bed > $LancetDir/$samplePair.NonLancet.called.bed && \\\n";
                $GetBedSH .= "cat $pairdir/$soft/$samplePair.Strelka2.indel.called.bed >> $LancetDir/$samplePair.NonLancet.called.bed && \\\n";
            }elsif( $soft eq "MuTect2" || $soft eq "MuSE" || $soft eq "SomaticSniper"){
            $GetBedSH .= "cat $pairdirShare/$soft/$samplePair.$soft.called.bed >> $LancetDir/$samplePair.NonLancet.called.bed && \\\n";
        }else{
                $GetBedSH .= "cat $pairdir/$soft/$samplePair.$soft.called.bed >> $LancetDir/$samplePair.NonLancet.called.bed && \\\n";
            }
    }
    $GetBedSH .= "sort -n -k 1.4,1.5 -k 2 -k 3 -u $LancetDir/$samplePair.NonLancet.called.bed | $bedtools intersect -a - -b $bed > $LancetDir/$samplePair.targetForLancet.bed && \\\n";
    my $LancetSH = $GetBedSH;
    $LancetSH .= "$Lancet --tumor $tbam --normal $nbam --ref $reffa --bed $LancetDir/$samplePair.targetForLancet.bed --num-thread 32 | $bgzip > $LancetDir/Lancet.$samplePair.raw.vcf.gz && \\\n";
    $LancetSH .= "zcat $LancetDir/Lancet.$samplePair.raw.vcf.gz|awk '\$1~/^#/ || \$7~/PASS/ {print}' | $bgzip > $LancetDir/$samplePair.Lancet.snv_shortindel.vcf.gz && \\\n";
    #$LancetSH .= &echoENDstring($shell06);
    if($local eq "T"){
    }else{
        open SH6, ">$shell06" || die $!;
        print SH6 &echoSTARTstring($shell06);
        print SH6 $LancetSH;
        print SH6 &echoENDstring($shell06);
        close SH6;
    }
    return $LancetSH;
}

##======================= calling somatic indel and sv by Svaba ==========================##
sub runSvaba{
    my ($pairdir,$samplePair,$nbam,$tbam,$shdir,$bed)=@_;
    my $SvabaDir = "$pairdir/Svaba";
    #unless(-e $SvabaDir){system("mkdir -p $SvabaDir")}
    my $shell08 = "$shdir/Svaba.$samplePair.sh";
    #open SH8, ">$shell08" || die $!;
    #print SH8 &echoSTARTstring($shell08);
    my $SvabaSH = "mkdir -p $SvabaDir && \\\n";
    $SvabaSH .= "$Svaba run -z -t $tbam -n $nbam -D $dbsnp -p 32 -a $SvabaDir/$samplePair -G $reffa && \\\n";
    $SvabaSH .= &getBed("$SvabaDir/$samplePair.svaba.somatic.indel.vcf.gz","$SvabaDir/$samplePair.Svaba.called.bed");
    #$SvabaSH .= &echoENDstring($shell08);
    if($local eq "T"){
    }else{
        open SH8, ">$shell08" || die $!;
        print SH8 &echoSTARTstring($shell08);
        print SH8 $SvabaSH;
        print SH8 &echoENDstring($shell08);
        close SH8;
    }
    return($SvabaSH);
}

##========================== sample QC for contaminations  ==========================##
## Determining whether two BAM contain reads sequenced from the same sample or patient by counting genotype matches at common SNPs.
sub somaticQC{
    my ($pairdir,$samplePair,$nbam,$tbam,$shdir)=@_;
    my $QCDir = "$pairdir/QC";
    unless(-e $QCDir){system("mkdir -p $QCDir")}
    system("mkdir -p $QCDir/tmp");
    my $shell10 = "$shdir/somaticQC.$samplePair.sh";
    open SH10, ">$shell10" || die $!;
    print SH10 &echoSTARTstring($shell10);
    #print SH10 "export PATH=$pythonBin:\$PATH && \\\n";
    print SH10 "export LD_LIBRARY_PATH=/opt/bioinfo/tools/samtools/xz/v5.2.3/lib:/opt/python/lib:\$LD_LIBRARY_PATH && \\\n";
    print SH10 "export PYTHONPATH=$bammathPATH/pythonLIB/lib/python2.7/site-packages && \\\n";
    print SH10 "python $bammathPATH/bam-matcher.py --bam1 $nbam --bam2 $tbam --output $QCDir/bam_matcher.report.txt --vcf $bammathPATH/1kg.exome.highAF.merge.hg38.vcf.gz --caller gatk --dp-threshold 15 --gatk-mem-gb 4 --gatk-nt 1 --reference $reffa --do-not-cache --experimental --scratch-dir $QCDir/tmp && \\\n";
    print SH10 "$java8 -Xms4g -Xmx4g -Djava.io.tmpdir=$QCDir/tmp -XX:MaxMetaspaceSize=128m -XX:-UseGCOverheadLimit -XX:ParallelGCThreads=1 -jar $GATK3 -T ContEst -R $reffa -I:eval $tbam -I:genotype $nbam --popfile $hapmapAF -L $hapmapAF -isr INTERSECTION -o $QCDir/ContEst.report.txt && \\\n";
    print SH10 "rm -rf $QCDir/tmp && \\\n";
    print SH10 &echoENDstring($shell10);
    close SH10;
}

## vcf to bed, for merge
sub getBed{
    my ($vcfFile,$bedFile) = @_;
    my $getbed = "$Vt decompose -s $vcfFile | $bcftools norm -f $reffa| $bcftools query -f '\%CHROM\\t\%POS\\t\%END\\t%REF\\t%ALT\\n'|awk '{if (length(\$4)==1 && length(\$5)==1) print \$1\"\\t\"\$2\"\\t\"\$3;else print \$1\"\\t\"\$2-10\"\\t\"\$3+10 }' >$bedFile && \\\n";
    return $getbed;
}

## echo start sign for monitor
sub echoSTARTstring{
    my $sh=shift;
    my $ostr="#!/bin/bash\n";
    $ostr.="set -euo pipefail\n";
    $ostr.="echo ==========start at : `date` on `hostname` ==========\n";
    return $ostr;
}

## echo end sign for monitor
sub echoENDstring{
    my $sh=shift;
    my $ostr="echo ==========end at : `date` on `hostname` ========== && \\\n";
    $ostr.="echo Still_waters_run_deep 1>&2\n";
    return $ostr;
}