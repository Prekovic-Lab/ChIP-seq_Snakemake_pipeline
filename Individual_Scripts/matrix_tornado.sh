computeMatrix reference-point --referencePoint center \
        -R file1.bed \ #(file2.bed, ..) 
        -S bigwig_cond1.bw bigwig_cond2.bw \
        -b 10000 -a 10000 \ # distance before and after reference point (tss/peak point)
        --skipZeros -o output_matrix_name.gz \
        --outFileNameMatrix output_tab_format.tab \ # This matrix can easily be loaded into R or other programs
        --outFileSortedRegions genes_output.bed \ # File name in which the regions are saved after skiping zeros or min/max threshold values.
        -p 15
