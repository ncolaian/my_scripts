#!/bin/bash
## Script that uses the metassemble.py wrapper to metassemble assembly A and B, using
## assembly A as the primary assembly, and a 2kb mate-pair simulated library.

### Generate configuration file that will be input for metassemble.py
name=$1
path=/lustre/scr/n/c/ncolaian/diginorm_assem

mkdir ${path}/${name}/MergeMetassemble

echo -e "\
############################################\n\
###   Metassemble A and B configuration file\n\
############################################\n\
[global]\n\
\n\
bowtie2_threads=4\n\
bowtie2_read1=${path}/${name}/split_pe_aa.1.new,${path}/${name}/split_pe_ab.1.new\n\
bowtie2_read2=${path}/${name}/split_pe_aa.2.new,${path}/${name}/split_pe_ab.2.new\n\
bowtie2_maxins=3000\n\
bowtie2_minins=1000\n\
\n\

\n\
mateAn_A=1300\n\
mateAn_B=2300\n\
\n\
[1]\n
\n\
fasta=${path}/${name}/metavelvet/meta-velvetg.contigs.fa\n\
ID=A\n\
\n\
[2]\n\
\n\
fasta=${path}/${name}/ab_metav/meta-velvetg.contigs.fa\n\
ID=B\n\
" > ${path}/${name}/MergeMetassemble/B.A.metassemble.config

### Run metassemble

bsub -M 200 -q bigmem -R "span[hosts=1]" -n 4 "python /proj/cdjones_lab/ncolaian/apps/Metassembler/src/metassemble.py --conf ${path}/${name}/MergeMetassemble/B.A.metassemble.config --outd ${path}/${name}/MergeMetassemble"  

