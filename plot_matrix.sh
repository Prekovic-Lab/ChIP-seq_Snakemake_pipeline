#!/bin/bash

plotHeatmap -m matrix_from_last_command.gz -o FileName.pdf \
        --colorMap Reds -min 0 -max 20 \
        --whatToShow "heatmap and colorbar"
