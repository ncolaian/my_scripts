# my_scripts
This is a collection of some of my scripts that I have made. 

metagenome_scripts: 
 Contains scripts that were created for projects related to analyzing metagenome data.
  - analyzing_28S_data.pl => This program is a beginning script to view taxanomic data
  - analyzing_new_metagenome_tax_info.pl => This program creates a csv text file that can be passed to metagenome_taxa_graphing.R to create graphs analyzing the Eukaryotic species in a metagenome sample
  - get_tax_info.pl => will retrieve the taxonomy information of genomes from the 97_otu_taxonomy.txt file. ( Path on killdevil is /proj/dangl_lab/data/databases/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt )
  - kog_db_specific.pl => will trim a full kog file into a kog file that only contains the kogs found in the database
  - updated_rRNA_pred_main => original program made by Scott Yourstone. This program will run rRNA prediction to find all the 28S sequences in a metagenomic assembly file.
 
phylogeny_scripts: 
 Contains scripts that were created for projects related to phylogeny creation.
  - bsub_blast.sh => This is a simple bash script that will run blast on all the assembled drosophila data
  - bsub_get_spades_sequences.sh => will get the portion of the sequence that came up as a hit in the blast file ( good for single hits )
  - get_comparitive_do_pos.pl => used to get the drop out positions from a list of genes from different genomes. This is to ensure that the portion of sequence being pulled out exists in all the genomes. 
  - get_dropout_areas.pl => will take a file that has the amount of hits each base has and determine if it's a dropout or not. It will then record the dropout positions within the sequence. The dropouts are calculated by any base that has less that 3x the STD of reads map back to it.
  - remove_dropouts.pl => once you have the positions of all the dropouts, you can use this program to create concatenated genes w/o any dropouts
  
r_scripts: 
 Contains any R scripts that I have created. Mostly scripts to make ggplot graphs.
  - metagenome_taxa_graphing.R => Creates graphs that display eukaryotic taxa information.
  
utility_scripts: 
 Contains some utility scripts that were created to do jobs that could apply to multiple different projects
  - combine_fastq_files.pl => will create interleaved fastq files from separated paired end files.
  - make_pe_tags_identical.pl => ensures the tags on paired end tags are identical. Important for some downstream analysis
  - fasta_file_separator.pl => creates individual fasta files from a file that contains many entries
  - reverse_dna_file.pl => the sequences will be transformed to the opposite strand. So it will be reversed and bases will be switched.
