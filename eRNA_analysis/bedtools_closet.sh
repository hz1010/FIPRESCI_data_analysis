##bulk ATAC-seq data analysis is following guidelines in https://informatics.fas.harvard.edu/atac-seq-guidelines.html. Here we have already get narrow_peaks and caculating their shortest distance to TSS

#!/bin/sh
cd ~/distal_peak
bedtools closest -a ~/Hela_peaks.narrowPeak -b ~/regions/tss_sorted.bed -d > Hela_dis2k.bed

