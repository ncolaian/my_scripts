#! usr/bin/bash

# Run metavelvet
name=$1
#file_directory=$2
#output_path=$3

#files to run metavelvet on
#make sure to change file name based upon full or digital normalization

pefile="/lustre/scr/n/c/ncolaian/diginorm_assem/${name}_rz/${name}_ab"
#sefile=$file_directory/$name/$name.raw.trimmo.fastq.se.keep
output_path="/lustre/scr/n/c/ncolaian/diginorm_assem/${name}_rz/ab_metav"
mkdir $output_path


# Run velveth

export OMP_NUM_THREADS=5

#this is for data with se files
#bsub -J "${name}1" -q bigmem -M 1000 -n 6 -R "span[hosts=1]" "velveth ${output_path} 31 -fmtAuto -shortPaired ${pefile} -short ${sefile}"

#for data without paired files
bsub -J "${name}1" -q bigmem -M 1000 -n 6 -R "span[hosts=1]" "velveth ${output_path} 31 -fmtAuto -shortPaired ${pefile}" 

bsub -J "${name}2" -w "${name}1" -q bigmem -M 300 -R "span[hosts=1]" "velvetg ${output_path} -exp_cov auto -ins_length 260"

bsub -J "${name}3" -w "${name}2" -q week -M 100 -R "span[hosts=1]" "meta-velvetg ${output_path} -ins_length 260 | tee ${output_path}logfile.txt"
