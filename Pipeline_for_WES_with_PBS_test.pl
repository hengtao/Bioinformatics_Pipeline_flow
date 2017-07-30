#!/usr/bin/perl -w

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
=head usage()

perl Pipeline_for_LungCancer_noArgument.pl  -sampleInfo sampleInfo
sampleInfo file contains sample names, one couple per line,such as 0116D01M200001A01_2101	0116D01M200001C01_2101

-h: display this help

eg: perl Pipeline_for_LungCancer_noArgument.pl -sampleInfo sampleInfo.txt
=cut



use Getopt::Long;
use File::Path;


GetOptions("h" => \$help,
           "sampleInfo=s" => \$sampleInfo,
);

if ($help ) {
	usage();
	exit(1);
}

my $outputPath="/mnt/gennlife/Ucloud_test/test_results";
my $HOME="/home/bio";
my $logPath="$outputPath/0-log";
my $toolsPath="/opt/mnt/bio_tools";
my $metadataPath="/mnt/gennlife/data/Bio_data";
my $Fq_position="/mnt/gennlife/copy_to_ucloud";
my $pbsPath="/mnt/gennlife/Ucloud_test/test_results/PBS";



my $commands='';
open(SAMINFO,"<$sampleInfo") or die "$!\n";

while (<SAMINFO>)
{
    chomp;
    my ${line}=$_;
    next if(${line} =~ /^\s+$/);
    next if(${line} =~ /^\#/);
       
    
    #$commands .= "perl /opt/mnt/cephfs/bio/Pipeline_for_WES.pl ${line}\n";
    my $normal = (split("\t",$line))[0]; 
    my $tumor = (split("\t",$line))[1];   
    print "$tumor\n";
    if(! -e "$pbsPath/${tumor}")
      {
          mkpath("$pbsPath/${tumor}",0644);
          if($@)
          {
                          print "Make path $pbsPath/${tumor} failed\n";
                          exit(1);
                }
      }

    $commands .= "perl /mnt/gennlife/Pipeline_for_WES.pl $normal & \n";
    $commands .= "perl /mnt/gennlife/Pipeline_for_WES.pl $tumor & \n";
    $commands .= "wait \n";
    $commands .= "java -jar /opt/mnt/bio_tools/GenomeAnalysisTK/GenomeAnalysisTK.jar -T MuTect2 -R $metadataPath/hg19/hg19.fa -I:tumor $outputPath/4-GATK/$tumor.marked.realn.recal.bam -I:normal $outputPath/4-GATK/$normal.marked.realn.recal.bam --dbsnp $metadataPath/Mutect/dbsnp_138.hg19.mutect.vcf --cosmic $metadataPath/gatk/hg19.convert_cosmic_v54_120711.vcf -o $outputPath/7-Mutect2/$tumor.mutect2.vcf \n";
    
    if($commands ne '')
    {
        open(PBS,'>',"$pbsPath/${tumor}/WES_${tumor}.pbs") or die "$!\tcan\'t write into this PBS file:$pbsPath/${tumor}/WES_${tumor}.pbs\n";
        print PBS << "EOF";
#!/bin/bash

#PBS -N WES_${tumor}
#PBS -l nodes=1:ppn=4
#PBS -l walltime=100:00:00
#PBS -o $pbsPath/${tumor}/WES_${tumor}.log
#PBS -e $pbsPath/${tumor}/WES_${tumor}.err
#PBS -m abe
#PBS -M wanghengtao\@gennlife.com


date
$commands
date
EOF
;
        close PBS;
        my $qsub= `qsub $pbsPath/${tumor}/WES_${tumor}.pbs`;
        #my $qsub = '';
        if ($qsub ne '')
        {
            print "Task WES_${tumor} was submitted in $qsub\n";
        }
        my $taskNum=`qstat -u bio |grep WES|wc -l`;
        while($taskNum>=80)
        {
            sleep (100);
            $taskNum=`qstat -u bio |grep WES|wc -l`; 
        }
    }
    $commands='';
}



sub usage
{
        die `pod2text $0`;
}








