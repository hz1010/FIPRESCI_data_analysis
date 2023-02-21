##count RNA reads for each distal peaks

#!/bin/bash
bedtools intersect -a ~/PolyT_unique_RNA.bed -b  enhancer_tss2k.bed -c > Overlap_PolyT.bed

