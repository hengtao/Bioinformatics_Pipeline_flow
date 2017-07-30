#!/usr/bin/perl 

#########################################
#Copyright (c) Gennlife
#Author: WANG Hengtao
#Email: wanghengtao@gennlife.com
#Version: 2.0
#Date: 2016.08.05
#Motified: 2016.08.08
#Description: 
########################################

=pod
=head1 Usage

perl Pipeline_for_LungCancer_noArgument.pl  the_samplename

-h: display this help

eg: perl Pipeline_for_LungCancer_noArgument.pl 0116D01M200001C01_2101
=cut

#use strict;
use Getopt::Long;
use File::Path;
#use PerlIO::gzip;

GetOptions("h" => \$help);

if (@ARGV<1 || $help ) {
	usage();
	exit(1);
}


my ($string1,$string2,$ID,$SM,$PU,$LB); 
$SM=$ID=$ARGV[0];
my $PL="Illumina";

my $outputPath="/mnt/gennlife/Ucloud_test/test_results";
my $HOME="/home/bio";
my $logPath="$outputPath/0-LOG";
my $toolsPath="/opt/mnt/bio_tools";
my $metadataPath="/mnt/gennlife/data/Bio_data";
my $Fq_position="/mnt/gennlife/copy_to_ucloud";



if(!-e "$logPath"){
	mkpath("$logPath",0644);
	if($@){
		print "Make path $logPath failed\n";
		exit(1);
	}
}

if(!-e "$outputPath/0-FASTQC"){
	mkpath("$outputPath/0-FASTQC",0644);
	if($@){
		print "Make path $outputPath/0-FASTQC failed\n";
		exit(1);
	}
}


if(!-e "$outputPath/1-SPLIT/$SM"){
	mkpath("$outputPath/1-SPLIT/$SM",0644);
	if($@){
		print "Make path $outputPath/1-SPLIT/$SM failed\n";
		exit(1);
	}
}

if(!-e "$outputPath/2-MAPPING/"){
	mkpath("$outputPath/2-MAPPING/",0644);
	if($@){
		print "Make path $outputPath/2-MAPPING/ failed\n";
		exit(1);
	}
}


if(!-e "$outputPath/3-PICARD/"){
	mkpath("$outputPath/3-PICARD/",0644);
	if($@){
		print "Make path $outputPath/3-PICARD/ failed\n";
		exit(1);
	}
}


if(!-e "$outputPath/4-GATK/"){
	mkpath("$outputPath/4-GATK/",0644);
	if($@){
		print "Make path $outputPath/4-GATK/ failed\n";
		exit(1);
	}
}

if(!-e "$outputPath/5-GERMLINE/"){
	mkpath("$outputPath/5-GERMLINE/",0644);
	if($@){
		print "Make path $outputPath/5-GERMLIN/ failed\n";
		exit(1);
	}
}


if(!-e "$outputPath/6-HsMetrics/"){
	mkpath("$outputPath/6-HsMetrics/",0644);
	if($@){
		print "Make path $outputPath/6-HsMetrics/ failed\n";
		exit(1);
	}
}

if(!-e "$outputPath/7-Mutect2/"){
	mkpath("$outputPath/7-Mutect2/",0644);
	if($@){
		print "Make path $outputPath/7-Mutect2/ failed\n";
		exit(1);
	}
}



my $FASTQC = "$toolsPath/FastQC/fastqc";
my $BWA = "$toolsPath/bwa-0.7.12/bwa";
my $SAMTOOLS = "$toolsPath/samtools-1.3/samtools";
my $PICARDJAR="$toolsPath/picard-tools-1.141/picard.jar";
my $GATKJAR="$toolsPath/GenomeAnalysisTK/GenomeAnalysisTK.jar";



my $hgRef="$metadataPath/hg19";
my $known1000G_indels="$metadataPath/gatk/1000G_phase1.indels.hg19.sites.vcf";
my $GoldStandard_indels="$metadataPath/gatk/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf";
my $dbSNP="$metadataPath/gatk/dbsnp_138.hg19.vcf";
my ${targets_interval_list}="$metadataPath/S07504514.interval_list";


my $file1=$ARGV[0]."_R1.fq.gz";
my $file2=$ARGV[0]."_R2.fq.gz";

if (-e "$Fq_position/$file1") 
{ 
    system "fastqc -o $outputPath/0-FASTQC  $Fq_position/${ARGV[0]}_R1.fq.gz";
    system "echo fastqc finished for $ARGV[0]  >>  $ARGV[0]\.log";
    open (IN1, "gzip -dc $Fq_position/$ARGV[0]_R1.fq.gz  | ") || die "$!\tthis file not exists!\n"; 
}else
{
    system "fastqc -o $outputPath/0-FASTQC  $Fq_position/$ARGV[0]/${ARGV[0]}_R1.fq";
    system "echo fastqc finished for $ARGV[0]  >>  $ARGV[0]\.log";
    open (IN1, "$Fq_position/$ARGV[0]/$ARGV[0]_R1.fq") || die "$!\tthis file not exists!!\n";
}

if (-e "$Fq_position/$file2") 
{ 
    system "fastqc -o $outputPath/0-FASTQC  $Fq_position/${ARGV[0]}_R2.fq.gz";
    system "echo fastqc finished for $ARGV[0]  >>  $ARGV[0]\.log";
    open (IN2, "gzip -dc $Fq_position/$ARGV[0]_R2.fq.gz  | ") || die "$!\tthis file not exists!!!\n"; 
}else 
{
    system "fastqc -o $outputPath/0-FASTQC  $Fq_position/$ARGV[0]/${ARGV[0]}_R2.fq";
    system "echo fastqc finished for $ARGV[0]  >>  $ARGV[0]\.log";
    open (IN2, "$Fq_position/$ARGV[0]/$ARGV[0]_R2.fq") || die "$!\tthis file not exists!!!!\n";
}


my $string="";



my $line1=<IN1>; my $line2=<IN2>;
chomp $line1; chomp $line2;
my @array1=split(/:/,$line1);
my @array2=split(/:/,$line2);

$DI=$array1[2].".".$array1[3];
$LB=$array1[2];
$PU=$DI.".".$SM;





my $new1=$DI."1"; my $new2=$DI."2";
open ($new1, ">$outputPath/1-SPLIT/$SM/${DI}.1.fq") || die "$!can not write!\n";
open ($new2, ">$outputPath/1-SPLIT/$SM/${DI}.2.fq") || die "$!can not write!!\n";

print $new1 "$line1\n";
print $new2 "$line2\n";

my $j=1;

while (defined $line1)
{
    $line1=<IN1>; 
    $line2=<IN2>;
    if (! defined $line1)
    {
       last;
    }
    chomp ($line1,$line2);
    
    if($line1=~/^\@/)
    {
        my $DI_old=$DI;
        my $LB_old=$LB;
        my $PU_old=$PU;
        
        @array1=split(/:/,$line1);
        @array2=split(/:/,$line2);

        $DI=$array1[2].".".$array1[3];
        $LB=$array1[2];
        $PU=$DI.".".$SM;
        if ($DI eq $DI_old)
        {
            $new1=$DI."1";
            $new2=$DI."2";
            
            print $new1 "$line1\n";
            print $new2 "$line2\n";
        }
        else 
        {
            
            close $new1;
	    close $new2;
            
            system "bwa mem -t 2 -R '\@RG\\tID:$DI_old\\tSM:$SM\\tPU:$PU_old\\tLB:$LB_old\\tPL:$PL' -M $hgRef/hg19.fa $outputPath/1-SPLIT/$SM/${DI_old}.1.fq   $outputPath/1-SPLIT/$SM/${DI_old}.2.fq | samtools view -bS -o $outputPath/2-MAPPING/${SM}.${DI_old}.bam -";
            #system "samtools view -bS -o $outputPath/2-MAPPING/${SM}.${DI_old}.bam $outputPath/2-MAPPING/${SM}.${DI_old}.sam >> $logPath/${SM}.${DI_old}.bwa.log";
            system "echo \"bwa finished for $SM $DI_old\" >>  $logPath/$SM\.log";
            system "echo `date` >> $logPath/$SM\.log";
            
            system "samtools sort  $outputPath/2-MAPPING/${SM}.${DI_old}.bam -T $outputPath/2-MAPPING/${SM}.${DI_old}.${j} -o $outputPath/2-MAPPING/${SM}.${DI_old}.sorted.bam";
            system "echo \"sort finished for $SM $DI_old\" >>  $logPath/$SM\.log";
            system "echo `date` >> $logPath/$SM\.log";
            			      
	    $new1=$DI."1";
            $new2=$DI."2";
            
            open ($new1, ">$outputPath/1-SPLIT/$SM/${DI}.1.fq") || die "$!\n";
            open ($new2, ">$outputPath/1-SPLIT/$SM/${DI}.2.fq") || die "$!\n";
            
            print $new1 "$line1\n";
            print $new2 "$line2\n";
            
            $string.=" I="."$outputPath/2-MAPPING/${SM}.${DI_old}.sorted.bam ";
            $j++;
         }
     }
     else
     {
        print $new1 "$line1\n";
        print $new2 "$line2\n";
     }
}

close $new1; close $new2; close IN1; close IN2;
 

system "bwa mem -t 2 -R '\@RG\\tID:$DI\\tSM:$SM\\tPU:$PU\\tLB:$LB\\tPL:$PL' -M $hgRef/hg19.fa $outputPath/1-SPLIT/$SM/${DI}.1.fq $outputPath/1-SPLIT/$SM/${DI}.2.fq |samtools view -bS -o $outputPath/2-MAPPING/${SM}.${DI}.bam -";
system "echo \"bwa finished for $SM $ID\" >>  $logPath/$SM\.log";
system "echo `date` >> $logPath/$SM\.log";
            
system "samtools sort  $outputPath/2-MAPPING/${SM}.${DI}.bam -T $outputPath/2-MAPPING/${SM}.${DI}.${j} -o $outputPath/2-MAPPING/${SM}.${DI}.sorted.bam";
system "echo \"sort finished for $SM $DI\" >>  $logPath/$SM\.log";
system "echo `date` >> $logPath/$SM\.log";

$string.=" I="."$outputPath/2-MAPPING/${SM}.${DI}.sorted.bam ";

system "echo  \"the $SM was split into $j files\" >> $logPath/$SM\.log";            
                   

#system "samtools merge -nufr -b $outputPath/sorted.bamlist.txt --threads 4 -O $outputPath/2-MAPPING/$SM.bam >> $logPath/$SM.Samtools_merge\.log";
system "java -jar $PICARDJAR MergeSamFiles $string OUTPUT=$outputPath/2-MAPPING/$SM.bam";
system "echo \"merge finished for $SM\" >>  $logPath/$SM\.log";
system "echo `date` >> $logPath/$SM\.log";



system "java -jar $PICARDJAR MarkDuplicates VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=True INPUT=$outputPath/2-MAPPING/$SM.bam OUTPUT=$outputPath/3-PICARD/$SM\.marked.bam METRICS_FILE=$outputPath/3-PICARD/$SM.Mkdup.metrics";
system "echo \"markduplicate finished for $SM\" >>  $logPath/$SM\.log";
system "echo `date` >> $logPath/$SM.log";

system "samtools index $outputPath/3-PICARD/$SM\.marked.bam";
system "echo \"samtools index finished for $SM\" >>  $logPath/$SM\.log";
system "echo `date` >> $logPath/$SM.log";



    
system "java -Xmx3g -jar $GATKJAR  -T RealignerTargetCreator -U ALLOW_N_CIGAR_READS -R $hgRef/hg19.fa -I $outputPath/3-PICARD/$SM.marked.bam -o $outputPath/4-GATK/$SM.marked.intervals -known $known1000G_indels -known $GoldStandard_indels --interval_padding 150 ";
system "echo \"GATK Realinger finished  for $SM\" >>  $logPath/$SM\.log";
system "echo `date` >> $logPath/$SM.log";

system "java -Xmx3g -jar $GATKJAR -T IndelRealigner -R $hgRef/hg19.fa -I $outputPath/3-PICARD/$SM.marked.bam -targetIntervals $outputPath/4-GATK/$SM.marked.intervals -o $outputPath/4-GATK/$SM.marked.realn.bam -known $known1000G_indels -known $GoldStandard_indels --interval_padding 150 -LOD 0.5";
system "echo \"GATK Indel Realinger finished  for $SM\" >>  $logPath/$SM\.log";
system "echo `date` >> $logPath/$SM\.log";

system "java -Xmx4g -jar $GATKJAR -nct 2 -T BaseRecalibrator -R $hgRef/hg19.fa -I $outputPath/4-GATK/$SM.marked.realn.bam -o $outputPath/4-GATK/$SM.marked.realn.recal -knownSites $known1000G_indels -knownSites $GoldStandard_indels -knownSites $dbSNP --interval_padding 150";
system "echo \"GATK BaseRecalibrator finished  for $SM\" >>  $logPath/$SM\.log";
system "echo `date` >> $logPath/$SM\.log";

system "java -Xmx3g -jar $GATKJAR  -T PrintReads -R $hgRef/hg19.fa -I $outputPath/4-GATK/$SM.marked.realn.bam -o $outputPath/4-GATK/$SM.marked.realn.recal.bam --BQSR $outputPath/4-GATK/$SM.marked.realn.recal --interval_padding 150";
system "echo \"GATK PrintReads finished  for $SM\" >>  $logPath/$SM\.log";
system "echo `date` >> $logPath/$SM\.log";

system "java -Xmx3g -jar $GATKJAR -T HaplotypeCaller -R $hgRef/hg19.fa -L ${targets_interval_list} --dbsnp ${dbSNP} --emitRefConfidence GVCF -I $outputPath/4-GATK/$SM.marked.realn.recal.bam -o $outputPath/5-GERMLINE/$SM.marked.realn.recal.bam.g.vcf --interval_padding 150 PCR=\"CONSERVATIVE\"";
system "echo \'GATK HaplotypeCaller finished  for $SM\' >>  $logPath/$SM\.log";
system "echo `date` >> $logPath/$SM\.log";

system "java -jar $PICARDJAR CalculateHsMetrics I=$outputPath/4-GATK/$SM.marked.realn.recal.bam O=$outputPath/6-HsMetrics/$SM.marked.realn.recal.bam.metrics BAIT_INTERVALS=${targets_interval_list} TARGET_INTERVALS=${targets_interval_list}";
system "echo \'GATK CalculateHsMetrics  finished  for $SM\' >>  $logPath/$SM\.log";
system "echo `date` >> $logPath/$SM\.log";



sub usage
{
        die `pod2text $0`;
}


### THis is the command line for this script###
#perl Pipeline_for_LungCancer.pl --sampleFile1 SampleFile_1.fq.gz --sampleFile2 SampleFile_2.fq.gz --toolsPath /opt/bio/bio_tools --metadataPath /opt/bio/bio_data --logPath /opt/mnt/vdb/log  --outputPath /opt/mnt/vdb
