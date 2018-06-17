#!/usr/bin/env python2
'''
Input: (1) Path to the directory storing all genomes (2) Path to the file containing the sqlite3 db.
This script gets a list of GOLD_IDs from the GOLD_FEATURES table. Then, for each GOLD_ID, it does the following:
It gets a list of GENBANK_IDs from the temp_acc table. For each GENBANK_ID, it finds the associated proteome file
and maintains a running count of the number of occurrences of IVYWREL and the total number of amino acids. It calculates
the fraction of IVYWREL amino acids and updates the GOLD_FEATURES table with this value.
Run as follows: ./calc_f_ivywrel_publish.py /path/to/genomes/directory /path/to/sqlite3/db/file
'''
import sys
import Bio
import numpy as np
from Bio import SeqIO
import os
from gzip import GzipFile
import sqlite3
import re


if __name__ == '__main__': 

    all_genomes_path = sys.argv[1]  # Path to directory storing all genomes. 
    db = sys.argv[2]    # Path to sqlite3 db. 
    conn = sqlite3.connect(db) 
    c = conn.cursor()

    # Get a list of GOLD_IDs from the db that have OGTs.
    query = "SELECT GOLD_ID FROM GOLD_FEATURES WHERE OGT IS NOT NULL;"
    c.execute(query)
    goldIDs = [i[0] for i in c.fetchall()]

    for goldID in goldIDs:
        
        ivywrel = re.compile(r'[IVYWREL]')   # Only IVYWREL amino acids.
        all_aa = re.compile(r'[GAPVLIMFYWSTCNQKHRDE]')    # All amino acids
        ivywrel_count = 0 # Stores the number of AAs that are belong to the IVYWREL set.
        total_length = 0 # Stores the total number of AAs. 

        # Get all GENBANK_IDs corresponding to goldID.
        query = "SELECT DISTINCT GENBANK_ID FROM temp_acc WHERE GOLD_ID='{}';".format(goldID)
        c.execute(query)
        genbankIDs = [i[0] for i in c.fetchall()]
        for genbankID in genbankIDs:
            # Get the file handle for the proteome corresponding to the genbankID:
            genome_path = os.path.join(all_genomes_path, genbankID)
            if not os.path.exists(genome_path):
                sys.stderr.write("Path does not exist: {}. Continuing to next genbank ID.".format(genome_path))
                continue
            else:
                proteome_path = None
                for candidate_filename in os.listdir(genome_path):
            	    if candidate_filename.endswith("protein.faa.gz"):
            		    proteome_path = os.path.join(genome_path, candidate_filename)
            		    break
                if proteome_path is None:
                    sys.stderr.write("No proteome file for genbank ID {}\n".format(genbankID))
                    continue
                elif proteome_path.endswith("gz"):
                    handle = GzipFile(proteome_path)
                else:
                    handle = file(proteome_path, 'r')

                # If proteome file exists, start counting IVYWREL and total number of AAs.
                for line in handle:
                    if line[0] == '>':  # Skip the header line.
                        continue
                    else:
                        ivywrel_count += len(ivywrel.findall(line))
                        total_length += len(all_aa.findall(line))

        if total_length == 0: # This will happen only if none of the genbankIDs in the previous step had proteome files.
            sys.stderr.write("{} has no proteomes\n".format(goldID))
        else:
            # Update the db
            f_ivywrel = ivywrel_count / float(total_length)
            query = "UPDATE GOLD_FEATURES SET F_IVYWREL='{}' WHERE GOLD_ID='{}';".format(f_ivywrel, goldID)
            c.execute(query)
            conn.commit()

    conn.commit()
    conn.close()

