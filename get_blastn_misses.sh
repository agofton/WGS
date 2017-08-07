#!/bin/bash
#Requires blastn and usearch on the PATH

#BLAST parameters
################################################################################
input_dir="" #complete path
output_dir="" #complete path , output name with be input_name.txt
db_dir="" #complete path
db_name=""
num_threads="" #integer ,  defult = 1
evalue="0.00000000001"
max_target_seqs="1"

ids_dir=""

blast_fasta_out_dir=""
################################################################################
# step 1 perform BLAST sheach

      mkdir ${output_dir}

for file1 in i${input_dir}/*.fasta
  do
    blastn -task blastn -query ${file1} -db ${db_dir}/${db_name} -evalue ${evalue} -max_target_seqs ${max_target_seqs} -num_threads ${num_threads} -outfmt 6 > "${output_dir}/$(basename "$file1" .fasta)_blastout.txt" 
done

#the was the blastn script is formatted with > output instead of -out is important for the blastn_progress.sh script
################################################################################
# step 2, extract hits seq IDs

      mkdir ${ids_dir}

for file2 in ${output_dir}/.txt
  do
    cat ${file2} | awk '{print $1;}' > "${ids_dir}/$(basename "$file2" .txt)_seqids.txt"
done
################################################################################
# step 3 filtering hits from raw data file by seq ids

      mkdir ${blast_fasta_out_dir}

for file3 in ${input_dir}/.fasta
  do
    usearch9.2 -fastx_getseqs ${file3} -labels "${ids_dir}/$(basename "$file3" .fasta)_blastout_seqids.txt" -notmatched "${blast_fasta_out_dir}/$(basename "$file3" .fasta)_blastout.fasta"
done
