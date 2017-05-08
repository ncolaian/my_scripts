#This will take a file that has sample #'s and produce a combined fastq file containing all the samples. The path is the barseq files path

file=$1
out=$2
path="/netscr/ncolaian/barseq/fastq/CL21_forward"

bsub_cmd="bsub \""

while read num
do
    if (( $num > 10));
    then
	add=" gunzip -c ${path}/CL21_IT0${num}_L${num}_R1.fastq.gz >> ${out};"
    else
	add=" gunzip -c ${path}/CL21_IT00${num}_L${num}_R1.fastq.gz >> ${out};"
    fi
    bsub_cmd=$bsub_cmd$add
done < $file

bsub_cmd="${bsub_cmd}\""

eval $bsub_cmd


    
	
