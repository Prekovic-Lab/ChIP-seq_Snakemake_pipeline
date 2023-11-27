#!/bin/bash

# for broad peak calling use the following
# modify the elements inside the brackets with your files
macs3 callpeak -t "input_bam" -c "control_bam" -n "output_name" --outdir peaks -f BAM -g hs --broad --broad-cutoff 0.1

#macs3 callpeak: This is the main command. macs3 is a program for the analysis of ChIP-seq data, and callpeak is a subcommand of macs3 that is used to call peaks from alignment results.
#-t "input_bam": The -t flag specifies the treatment file(s). Here, it’s set to the value of the input_bam variable, which should be the path to a BAM file containing your ChIP-seq data.
#-c "control_bam": The -c flag specifies the control file(s). Here, it’s set to the value of the control_bam variable, which should be the path to a BAM file containing your control data.
#-n "output_name": The -n flag specifies the name string of the experiment. Here, it’s set to the value of the output_name variable.
#--outdir peaks: The --outdir flag specifies the directory where the output files will be saved. Here, it’s set to peaks. If the folder does not exist it will be generated.
#-f BAM: The -f flag specifies the format of the input files. Here, it’s set to BAM.
#-g hs: The -g flag specifies the effective genome size. Here, it’s set to hs, which stands for Homo sapiens (human).
#--broad: This flag tells macs3 to call broad regions. This is typically used for histone modifications.
#--broad-cutoff 0.1: The --broad-cutoff flag specifies the cutoff for the broad region. Here, it’s set to 0.1.

# for narrow peak calling use the following
# modify the elements inside the brackets with your files
macs3 callpeak -t "input_bam" -c "control_bam" -n "output_name" --outdir peaks -f BAM -g hs -B -q 0.01

#-B: This flag tells macs3 to output bedGraph files for the pileup of shifted tags for treatment and control.
#-q 0.01: The -q flag specifies the q-value cutoff to call significant regions. Here, it’s set to 0.01. This is kind of default.
