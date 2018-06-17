#!/usr/bin/env python2

'''The script goes to the table and finds all the Genbank IDs associated \
with each GOLD ID. Then for each Genbank ID, it goes to the genbank ID dir, opens the rna \
file and grabs all the 16S sequences. It combines all the 16S sequences for all the Genbank \
IDs and combines them into a file and passes it to SINA. The output of SINA is an alignment \
for each seq. That output is then used to calculate f_gc.
Input: (1) Path to genomes dir.
(2) Path to dir that will store the files with combined 16S sequences.
(3) Path to dir that will store the MSA files.
(4) Path to db. 
(5) gold ID
(6) Path to helix file.
'''

import sys
import Bio
from Bio import SeqIO
import numpy as np
import os
from os import path
from gzip import GzipFile
from subprocess import check_call
import sqlite3

all_genomes_path = sys.argv[1]
output_dest_fna = sys.argv[2]
output_dest_msa = sys.argv[3]
db = sys.argv[4]
goldID = sys.argv[5]
conn = sqlite3.connect(db)
c = conn.cursor()

Gp_id_tup = (goldID,) 
ALN_WIDTH=50000 # Length of the SILVA alignment.
helix_path = sys.argv[6] 

# Load helix struc into a numpy array    

helix=""
f = open(helix_path)
trash = f.readline()  # Remove header line
for line in f:
	helix+=line.strip().replace(" ","")  
assert len(helix) == ALN_WIDTH, "Length of helix was %s (expected 50,000)"%(len(helix))

helix = helix+'.'  # Add a character because SINA returns alignments that are 50,001 characters long.

helix_array = np.array([list(helix)]) 
h_set = (">", "<", "[", "]")
helix_mask = np.in1d(helix_array, h_set)

# Pull 16S sequences out of each "genome" (genbank ID) of one "organism" (Gp ID)
rRNA_16S_handle = None 
rna_16S_filename = Gp_id_tup[0]+".16S.fna"
rRNA_16S_path = os.path.join(output_dest_fna, rna_16S_filename) 
rRNA_16S_handle = file(rRNA_16S_path, 'w')

for item in c.execute("SELECT DISTINCT GENBANK_ID FROM temp_acc WHERE GOLD_ID='{}';".format(goldID)):
    genome_path = path.join(all_genomes_path, item[0])
    rna_path = None
    for candidate_filename in os.listdir(genome_path):
	    if candidate_filename.endswith("rna_from_genomic.fna.gz"):
		    rna_path = os.path.join(genome_path, candidate_filename)
		    break
    if rna_path is None:
        raise ValueError, "No rna file for genbank ID {}".format(genbank)
    
    handle = GzipFile(rna_path)

    for record in SeqIO.parse(handle, "fasta"):
	    if "16S ribosomal RNA" in record.description:
		    rRNA_16S_handle.write(record.format('fasta'))

rRNA_16S_handle.close()

# Run SINA
rRNA_16S_msa_filename = Gp_id_tup[0]+".16S.msa.fasta"
rRNA_16S_msa_path = os.path.join(output_dest_msa, rRNA_16S_msa_filename)
command = "sina -i {} --intype=fasta -o {} --outtype=fasta --ptdb /home/data/silva/SSURef_NR99_128_SILVA_07_09_16_opt.arb".format(rRNA_16S_path, rRNA_16S_msa_path) 
print command
check_call(command, shell=True, stdout=sys.stderr, stderr=sys.stderr)

# Open MSA from SINA.
handle = open(rRNA_16S_msa_path)
l =[]
for record in SeqIO.parse(handle, "fasta"):
    print record
    l.append(list(record))
#print l
seq_array = np.array(l) # This should now be an array in which each sequence makes up one \
# line and the chars of the seq are the elments of that line. This should have 50,001 chars.


gc_sum = (seq_array == 'C').sum(axis=0)+(seq_array == 'G').sum(axis=0)
gc_sum_helix = gc_sum[helix_mask].astype(np.float)
at_sum = (seq_array == 'A').sum(axis=0)+(seq_array == 'T').sum(axis=0)
at_sum_helix = at_sum[helix_mask].astype(np.float)
frac_gc = str(np.nanmean(gc_sum_helix / (at_sum_helix + gc_sum_helix)))

# print "UPDATE GC_RNA SET f_gc = {} WHERE GOLD_ID = '{}';".format(frac_gc, goldID)
c.execute("UPDATE GC_RNA SET f_gc = {} WHERE GOLD_ID = '{}';".format(frac_gc, goldID)
conn.commit()
conn.close()

 



    

    
