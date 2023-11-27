#!/bin/bash

bamCoverage --bam merged.bam -o merged.bw \
    --binSize 10
    --normalizeUsing RPGC
    --effectiveGenomeSize 2150570000
    --extendReads


#--bam merged.bam: This flag specifies the input BAM file to process.
#-o merged.bw: This flag specifies the output file name.
#--binSize 10: This flag sets the size of the bins for the genome coverage calculation. In this case, the bin size is set to 10 base pairs. Smaller bin size for higher resolution but it will slower.
#--normalizeUsing RPGC: This flag sets the method to normalize the number of reads per bin. RPGC stands for Reads Per Genome Coverage.
#--effectiveGenomeSize 2150570000: This flag sets the effective genome size, which is used for normalization. This is for mouse example. See more for the rest: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
#--extendReads: This flag allows the extension of reads to fragment size. It’s generally not recommended for spliced-read data, such as RNA-seq, as it would extend reads over skipped regions.
