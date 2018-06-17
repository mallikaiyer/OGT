#!/bin/bash

# Directory names passed as arguments should not end with a '/'
genbankID_file=${1} # Path to file containing list of distinct genbank IDs.
ten_genes_dir=${2} # Path to dir containing hmms for ten proteins.
hmmsearch_output_dir=${3} # Path to dir that will store the hmmsearch output files.
genomes_dir=${4} # Path to dir containing all genomes.
for i in `less ${genbankID_file}`;
do
	mkdir ${hmmsearch_output_dir}/${i}
	for j in `ls ${ten_genes_dir}`;
	do
		hmmsearch --cut_tc --tblout=${hmmsearch_output_dir}/${i}/${j##*/}.out ${ten_genes_dir}/${j##*/} ${genomes_dir}/${i}/*_protein.faa* 
	done
done

