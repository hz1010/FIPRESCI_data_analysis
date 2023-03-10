import os
import time
import pandas as pd
t1=time.time()

##read round1 index
index=pd.read_table("../round1.txt",header=None)
data_R1=index.copy()
data_R2=index.copy()

##read fastq for demultiplexing
data1=open("1_R1_S.fq","r")
data2=open("1_R2_S.fq","r")
length_index=len(index)

for k in range(length_index):
    data_R1[0][k]=open(index[0][k]+"1_R1.fq","w")
    data_R2[0][k]=open(index[0][k]+"1_R2.fq","w")
lines1=data1.readlines()
lines2=data2.readlines()
length=len(lines1)
i=1

while i < length-1:
    for j in range(length_index):
        if lines2[i][0:6]==index[0][j]:
            data_R1[0][j].write(lines1[i-1])
            data_R1[0][j].write(lines1[i])
            data_R1[0][j].write(lines1[i+1])
            data_R1[0][j].write(lines1[i+2])
            data_R2[0][j].write(lines2[i-1])
            data_R2[0][j].write(lines2[i])
            data_R2[0][j].write(lines2[i+1])
            data_R2[0][j].write(lines2[i+2])

    i=i+4
data1.close()
data2.close()
for k in range(length_index):
    data_R1[0][k].close()
    data_R2[0][k].close()

print("success")
t2=time.time()
print(t2-t1)

##if you fastq files is big, please split them into small fastq files. Or you can use seqkit (https://bioinf.shenwei.me/seqkit/)
