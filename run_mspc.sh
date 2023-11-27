#!/bin/bash

mspc -i file1.narrowPeak file2.narrowPeak file3.narrowPeak -r bio -w 1e-4 -s 1e-8 -a 0.05 -p config.json --excludeHeader

#-i: This flag is used to specify the input files. In this case, file1.narrowPeak, file2.narrowPeak, and file3.narrowPeak are the input files.
#-r: This flag is used to specify the type of replicates. The bio value indicates that the replicates are biological.
#-w: This flag is used to specify the weak threshold. The value 1e-4 is the threshold set for weak peaks in this command.
#-s: This flag is used to specify the strong threshold. The value 1e-8 is the threshold set for strong peaks in this command.
#-a: This flag is used to specify the alpha value, which is the maximum allowable error rate. The value 0.05 indicates that the maximum error rate is 5%.
#-p: This flag is used to specify the output file for the peak calling parameters. In this case, config.json is the output file.
#--excludeHeader: This flag indicates that the header should be excluded from the output.
