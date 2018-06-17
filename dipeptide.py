#!/usr/bin/env python2

"""
Script starts by selecting a gold_id from DIPEPTIDE_FREQUENCIES. It then creates a file containing 
all the protein records for that gold_id, from all genbank_id associated with it. It then reads each
record in the file and calculates the frequency of all possible dipeptides. The DIPEPTIDE_FREQUENCIES
table updated with this number.

Input: (1) Path to genomes dir.
(2) Path to dir which will store combined proteomes
(3) Path to db

Run as follows:
./dipeptide.py /path/to/genomes/dir/ /path/to/combined/proteomes/dir/ /path/to/db/

"""

import sys
import Bio
import numpy as np
from Bio import SeqIO
import os
from os import path
from gzip import GzipFile
import sqlite3

all_genomes_path = sys.argv[1] # "/home/iyer/research/genomes_ncbi/genomes"
output_dest = sys.argv[2] # "/home/iyer/research/dipeptide/analyses/combined_proteomes/"
db = sys.argv[3]
conn = sqlite3.connect(db)
c = conn.cursor()

query = "SELECT GOLD_ID FROM DIPEPTIDE_FREQUENCIES;"
c.execute(query)
results = c.fetchall()
for result in results:
    goldID = result[0]
    print "Working on GOLD_ID: ", goldID
    filename = goldID+".comb.faa"
    organism_proteome_path = os.path.join(output_dest, filename)
    organism_proteome_handle = open(organism_proteome_path, 'w')

    query = "SELECT DISTINCT GENBANK_ID FROM temp_acc WHERE GOLD_ID='{}';".format(goldID)
    for item in c.execute(query):
        candidate_path = path.join(all_genomes_path, item[0])
        if os.path.exists(candidate_path):
            genome_path = candidate_path
        else:
            print candidate_path, " does not exist."
            continue
        proteome_path = None
        for candidate_filename in os.listdir(genome_path):
            if candidate_filename.endswith(("protein.faa.gz", "protein.faa")):
                proteome_path = os.path.join(genome_path, candidate_filename)
                break
            else:
                continue

        if proteome_path is None:
            sys.stderr.write("No proteome file for genbank ID {}".format(item[0]))
        elif proteome_path.endswith('.gz'):
            handle = GzipFile(proteome_path)
        else:
            handle = file(proteome_path, 'r')

        for record in SeqIO.parse(handle, "fasta"):
            organism_proteome_handle.write(record.format('fasta'))

    if organism_proteome_handle.tell() == 0:
        sys.stderr.write("{} has no proteomes".format(goldID))
    else:
        organism_proteome_handle.close()
        aa_list = ['G', 'A', 'I', 'L', 'V', 'M', 'C', 'D', 'E', 'N', 'Q', 'K', 'R', 'H', 'Y', 'W', 'F', 'S', 'T', 'P']
        aa_count_dict = {}
        aa_freq_dict = {}
        for i in aa_list:
            for j in aa_list:
                aa_count_dict[i+j] = 0
                aa_freq_dict[i+j] = 0

        for record in SeqIO.parse(organism_proteome_path, "fasta"):
            for i in range(len(record.seq)-1):
                dipeptide = str(record.seq[i:i+2])
                if dipeptide in aa_count_dict.keys():
                    aa_count_dict[dipeptide] += 1
                else:
                    continue

    total = 0
    for dipep in aa_count_dict.keys():
        total += aa_count_dict[dipep]
    for dipep in aa_freq_dict.keys():
        fraction = float(aa_count_dict[dipep]) / float(total)
        aa_freq_dict[dipep] = fraction
        query = "UPDATE DIPEPTIDE_FREQUENCIES SET F_{} = {} WHERE GOLD_ID = '{}';".format(dipep, fraction, goldID)
        c.execute(query)
        conn.commit()

conn.commit()
conn.close()
