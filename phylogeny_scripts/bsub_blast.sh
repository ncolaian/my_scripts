#fasta_dir=/proj/cdjones_lab/projects/waddell/jeremy/data

#for fa in albost_Luzon_modified.fasta albostrigata_SriLanka_modified.fasta bilim_Guam_modified.fasta bilim_Oahu_modified.fasta hypo_Luzon_modified.fasta immigrans_modified.fasta kep_Brunei_modified.fasta kep_Sarawak_modified.fasta koh_Phillipines_modified.fasta koh_Sarawak_modified.fasta nas_Mombasa_modified.fasta nas_Mysore_modified.fasta neon_Mysore_modified.fasta nivei_modified.fasta siam_Cambodia_modified.fasta sulf_NGuinea_modified.fasta sulf_NIreland_modified.fasta taxonf_modified.fasta taxong_modified.fasta taxonj_modified.fasta
#do
#    bsub -o blsf.out blastn -db blast_db/pallidi_dhd -query $fasta_dir/$fa -out blast_results_from_ass/$fa.blast.txt -outfmt "6 qseqid sseqid qstart qend score length nident" -max_target_seqs 25
#done


spades_dir=/proj/cdjones_lab/projects/waddell/assembly/spades_assemblies
abyss_dir=/proj/cdjones_lab/projects/waddell/assembly/abyss_assemblies

for fa in albom_Ishigaki albost_Borneo albost_Indonesia albost_Luzon albostrigata_SriLanka bilim_Guam bilim_Oahu hypo_Guam hypo_Luzon kep_Sarawak koh_Phillipines nas_Mombasa nas_Mysore neohypo_NGuinea neon_Mysore nivei pallidi pula_Sarawak siam_Cambodia sulf_NGuinea sulf_NIreland taxonJ
do
    bsub -o spades.out blastn -db blast_db/pallidi_dhd -query $spades_dir/$fa/contigs.fasta -out b_res_spades/$fa.blast.txt -outfmt "6 qseqid sseqid qstart qend score length nident" -max_target_seqs 25
done

for na in albost_Cambodia immigrans kep_Brunei koh_Sarawak taxonF taxonG 
do
    bsub -o abyss.out blastn -db blast_db/pallidi_dhd -query $abyss_dir/$na.fasta -out b_res_spades/$na.blast.txt -outfmt "6 qseqid sseqid qstart qend score length nident" -max_target_seqs 25
done
