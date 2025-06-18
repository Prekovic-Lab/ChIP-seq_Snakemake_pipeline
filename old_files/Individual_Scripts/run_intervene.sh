#!/bin/bash

intervene venn -i cond1_ConsensusPeaks.bed cond2_ConsensusPeaks.bed --save-overlaps

#intervene venn: This is the main command. Intervene is the tool, and venn is the module within the tool that is used to compute Venn diagrams of up to 6 sets.
#-i cond1_ConsensusPeaks.bed cond2_ConsensusPeaks.bed: The -i flag is used to specify the input files. Modify the input files.
#--save-overlaps: This flag tells Intervene to save overlapping regions/names for all the combinations as bed/txt files.
