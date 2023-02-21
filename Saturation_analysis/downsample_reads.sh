##pre-calculation of downsample size based on cell number and reads number
size=(14150000 12735000 11320000 9905000 8490000 7075000 5660000 4245000 2830000 1415000)

##One index for example
sample=TCAGTG
cd ~/${sample}
for i in {0..9}
do
echo ${size[i]}
echo ~/${sample}
mkdir size_${i}
##########script for dsub
cat > ~/size_${i}/split_${sample}_${i}.sh<<EOF
#!/bin/bash
#PBS -N split_${sample}_${i}
#PBS -q core40
#PBS -l mem=20gb,walltime=99:00:00,nodes=1:ppn=1
#PBS -o split_${sample}_${i}.log
#PBS -e split_${sample}_${i}.err
#PBS -V
#HSCHED -s Project_name+Software_Name+Species
##HSCHED -s hschedd default -a 500, -b 1500

cd ~/${sample}
seqtk sample -s100 ${sample}_S1_L001_R1_001.fastq ${size[i]} > ~/size_${i}/${sample}_${i}_S1_L001_R1_001.fastq
seqtk sample -s100 ${sample}_S1_L001_R2_001.fastq ${size[i]} > ~/size_${i}/${sample}_${i}_S1_L001_R2_001.fastq
EOF
dsub ~/size_${i}/split_${sample}_${i}.sh
done

