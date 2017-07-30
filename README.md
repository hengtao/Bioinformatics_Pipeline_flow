## Bioinformatics_Pipeline_flow

**This program will contain all the scripts that I use for data analysis, and will continuous updating** 

**1. WES pipeline -- Pipeline_for_WES.pl**   
```
This script was writen with perl, invoke some bash commands in it.
From fq format file, to recal.bam file all the way.
```
**2. PBS script -- Pipeline_for_WES_with_PBS_test.pl**
```
This script will invoke WES pipeline, and also run mutect2 for one paired sample.
Usage:
perl Pipeline_for_WES_with_PBS_test.pl --sampleInfo sampleInfo.txt
```
**3. Configure file -- sampleInfo.txt**
```
This is the configure file which contain the sample id of all samples you cope with， one couple per line.such as:
0116D01M200001A01_2101  0116D01M200001C01_2101
```


