computeMatrix reference-point --referencePoint center \
        -R file1.bed file2.bed \        #file3.bed, .. usually, we go with intervene output (lost, gained, common), this is for the regions that will be plotted (y-axis)
        -S bigwig_cond1.bw bigwig_cond2.bw \        # bigwig with the reads data (x-axis) [One bw file is one tornado plot]
        -b 10000 -a 10000 \        # distance before and after reference point (TSS/peak point) bigger: takes more time
        --skipZeros -o output_matrix_name.gz \        # the name of the table that will be used to plot the tornados/heatmaps  
        --outFileNameMatrix output_tab_format.tab \        # This matrix can easily be loaded into R or other programs
        --outFileSortedRegions genes_output.bed \        # File name in which the regions are saved after skipping zeros or min/max threshold values.
        -p 10        # number of processors (windows laptops at work usually have 12)

# the default binning is 10 with taking the mean reads
# you can change the binning with -bs or --binSize
# above we have put --skipZeros which refers to whether regions with only scores of zero should be included or not. The default is to include them.
