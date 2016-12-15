spades_dir=/proj/cdjones_lab/projects/waddell/assembly/spades_assemblies
abyss_dir=/proj/cdjones_lab/projects/waddell/assembly/abyss_assemblies

for fa in albom_Ishigaki albost_Borneo albost_Indonesia albost_Luzon albostrigata_SriLanka bilim_Guam bilim_Oahu hypo_Guam hypo_Luzon kep_Sarawak koh_Phillipines nas_Mombasa nas_Mysore neohypo_NGuinea neon_Mysore nivei pallidi pula_Sarawak siam_Cambodia sulf_NGuinea sulf_NIreland taxonJ
do
    bsub -o seq.out "perl /nas02/home/n/c/ncolaian/scripts/pull_out_bres_portion.pl -btf b_res_spades/${fa}.blast.txt -af ${spades_dir}/${fa}/contigs.fasta -o spades_dhd_sequences/${fa}.dhd -name ${fa}"
done

for na in albost_Cambodia immigrans kep_Brunei koh_Sarawak taxonF taxonG 
do
    bsub -o seq.out "perl /nas02/home/n/c/ncolaian/scripts/pull_out_bres_portion.pl -btf b_res_spades/${na}.blast.txt -af ${abyss_dir}/${na}.fasta -o spades_dhd_sequences/${na}.dhd -name ${na}"
done
