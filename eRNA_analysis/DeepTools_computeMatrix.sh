##using distal peal center as reference point to calculate region V.S. Reads counts matrix

computeMatrix reference-point --referencePoint center -S ${Con[i]}_unique_R1.bigwig \
-R ../enhancer_tss2k.bed \
-o ${Con[i]}_unique_2k.mat.gz -a 2000 -b 2000 -p 12

