# my_scripts
This is a collection of some of my scripts that I have made. 

metagenome_scripts: 
 Contains scripts that were created for projects related to analyzing metagenome data.
 
 
phylogeny_scripts: 
 Contains scripts that were created for projects related to phylogeny creation.
  
r_scripts: 
 Contains any R scripts that I have created. Mostly scripts to make ggplot graphs.
  - metagenome_taxa_graphing.R => Creates graphs that display eukaryotic taxa information.
  
utility_scripts: 
 Contains some utility scripts that were created to do jobs that could apply to multiple different projects
  - combine_fastq_files.pl => will create interleaved fastq files from separated paired end files.
  - make_pe_tags_identical.pl => ensures the tags on paired end tags are identical. Important for some downstream analysis
  - fasta_file_separator.pl => creates individual fasta files from a file that contains many entries
  - reverse_dna_file.pl => the sequences will be transformed to the opposite strand. So it will be reversed and bases will be switched.
  - pre_mult_align_seq_flip.pl => this ensure all the sequences in a multiple alignment file are on the same strand. Important to do before running muscle because it doesn't account for this.

lib:
 Contains some basic modules that are used in some of my scripts (A lot of them need to be updated)
 - Master_aln.pm => This module needs to be updated. Currently it creates a scoring matrix between two sequences. Can then get the location and value of the max score


Barseq_scripts:
This directory contains all my scripts that I created for the analysis of the barseq data. This does not contain the R code I used for the statistics. That is located on my personal computer. This does contain all the data processing though to run the barseq software in many different ways.