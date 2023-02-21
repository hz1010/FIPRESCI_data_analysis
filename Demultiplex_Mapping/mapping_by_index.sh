##mapping by cellranger for each round1 index, you can replace index by reading index_files especially you have a bigger number of indexes (384). The following script shows 28 indexes.

index=(AAGTAT AGATGT ATCCGG CAGGAG CGCATA CTTACG GATCAA AATTGG AGCACG ATGAAG CATAGA CGGCGT CTTGAA GCCAGA ACAAGG AGGTTA ATTAGT CCACGC CGGTCC GAAATA GCCGTT ACCCAA AGTAAA CAACCG CCGATG CGTTAT GAAGGG GCGAAT)

for i in {0..27}
do
cd ~/Cellranger
mkdir ${index[i]}
cd ${index[i]}
cat > Cellranger_${index[i]}.sh<< EOF
#!/bin/bash
#PBS -N SuperGirl_${i}
#PBS -q core24
#PBS -l mem=31gb,walltime=99:00:00,nodes=1:ppn=12
#PBS -o CellRanger_${index[i]}.log
#PBS -e CellRanger_${index[i]}.err
#PBS -V
#HSCHED -s Project_name+Software_Name+Species
##HSCHED -s hschedd

cd ~/Cellranger/${index[i]}
cellranger count --id=${index[i]} \
--fastqs=~/Group_by_index_fq/${index[i]} \
--sample=${index[i]} \
--chemistry="fiveprime" \
--localmem=20 \
--localcores=12 \
--include-introns \
--transcriptome=~/refdata-cellranger-mm10-3.0.0
EOF
dsub Cellranger_${index[i]}.sh
done


