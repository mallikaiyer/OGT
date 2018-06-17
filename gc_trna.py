#!/usr/bin/env python2

'''
Input: 
(1) Path to directory containing all genomes.
(2) Path to db.

The script goes to the table and gets a list of Gold IDs. For each 
GOLD_ID, it then gets a list of GENBANK_IDs. From the rna file for each 
genbank ID, it grabs all tRNA sequences, and keeps a running count of 
number of G+C bases and number of A+T bases. Afer looping through all
genbank IDs, it calculates the fraction of GC bases for the gold ID, and 
writes it to the DB.

Run as follows: ./gc_trna_publish.py /path/to/genomes/dir/ /path/to/db/
''' 

import sys
import Bio
from Bio import SeqIO
import os
from os import path
from gzip import GzipFile
import sqlite3
import numpy as np

all_genomes_path = sys.argv[1] # Path to dir containing all genomes.
db = sys.argv[2] # Path to db.
conn = sqlite3.connect(db) 
c = conn.cursor()

if __name__ == '__main__':

    base_set = ("A", "T", "G", "C")
    gc_set = ("G", "C")

    query = "SELECT GOLD_ID FROM GOLD_FEATURES WHERE OGT IS NOT NULL;"
    c.execute(query)
    goldIDs = [i[0] for i in c.fetchall()]
    for goldID in goldIDs:
        print "Working on goldID: {}".format(goldID)
        base_count = 0
        gc_count = 0
        query = "SELECT DISTINCT GENBANK_ID FROM temp_acc WHERE GOLD_ID='{}';".format(goldID)
        c.execute(query)
        genbankIDs = [i[0] for i in c.fetchall()]
        for genbankID in genbankIDs:
            print "Working on genbankID: {}".format(genbankID)
            genome_path = path.join(all_genomes_path, genbankID)
            if not os.path.exists(genome_path):
                    sys.stderr.write('Path does not exist: {}. Continuing to next genbank ID.'.format(genome_path))
                    continue
            else:
                rna_path = None
                for candidate_filename in os.listdir(genome_path):
                    if candidate_filename.endswith("rna_from_genomic.fna.gz") or candidate_filename.endswith("rna_from_genomic.fna"):
                        rna_path = os.path.join(genome_path, candidate_filename)
                        break
                if rna_path is None:
                    sys.stderr.write("No rna file for genbank ID {}\n".format(genbankID))
                    continue
                elif rna_path.endswith('.gz'):
                    handle = GzipFile(rna_path)
                else:
                    handle = file(rna_path, 'r')

                for record in SeqIO.parse(handle, "fasta"):
                    if "tRNA" in record.description:
                        for i in record:
                            if i in base_set:
                                base_count += 1
                            if i in gc_set:
                                gc_count += 1
                            else:
                                continue    	
                
                handle.close()

        if base_count == 0:
            sys.stderr.write("{} has no tRNA sequences.\n".format(goldID))
            continue
        else:
            frac_gc = gc_count/float(base_count)
            query = "UPDATE GOLD_FEATURES SET tRNA_GC = {} WHERE GOLD_ID = '{}';".format(frac_gc, goldID)
            c.execute(query)
            conn.commit()

    conn.close()


 



    

    


