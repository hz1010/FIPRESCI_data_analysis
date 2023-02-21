##following guidelines in https://github.com/seanken/FivePrime/tree/master/PeakCallingPipeline/paraclu

MinVal=1
PARACLU_OUT=~/cluster_p
min_density_rise=2.4
min_pos_with_data=0
min_sum=36
CTSS=~/paraclu.ctss
bamfile=cluster_p.bam
echo "Make CTSS"
make_ctss3.sh $bamfile $CTSS
echo "Run Paraclu"
RunParaclu.sh $CTSS $MinVal $PARACLU_OUT $min_density_rise $min_pos_with_data $min_sum



