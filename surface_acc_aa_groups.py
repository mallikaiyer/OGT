#!/usr/bin/env python2
'''
Input: 
(1) Path to db
(2) Path to dir with modpipe output files.

This script pulls UNQ_ID and BEST_MODEL from the BEST_MODELS table. For each best model, a structure is created and 
passed to the Bio.PDB DSSP function. The output of the function is parsed to determine which residues have relative
surface accessibility >= 0.05 (exposed). This list of residues is then used to calculate the frac of exposed 
positively-charged, negatively-charged, polar and hydrophobic amino acids. The AA_SURFACE_ACC table is updated with 
these fractions.

Run as follows: ./surface_acc_aa_groups_publish.py /path/to/db/ /path/to/modpipe/output/dir/
'''

import sys
import sqlite3
import os
from Bio.PDB import DSSP
from Bio.PDB.PDBParser import PDBParser
from gzip import GzipFile
from collections import Counter

DB = sys.argv[1] # Path to db 
MODELS_DIR = sys.argv[2] # Path to dir with modpipe output files 
RSA_EXPOSED_THRESH = 0.05 # above this fraction of relative surface accessability, consider residue exposed
conn = sqlite3.connect(DB)
c = conn.cursor()

query = "select UNQ_ID, BEST_MODEL from BEST_MODELS;"
c.execute(query)
results = [i for i in c.fetchall()]
for result in results: 
    unq_id    = result[0]
    bestmodel = result[1]
    print "This is the bestmodel {} for unq_id {}".format(bestmodel, unq_id)

    filename = os.path.join(MODELS_DIR, unq_id[:3], unq_id, "models", "{}.pdb.gz".format(bestmodel))
    if os.path.exists(filename):
        pdb_path = filename
    elif os.path.exists(filename[:-3]): # if the file is not gzipped
        pdb_path = filename[:-3]
    else:
        sys.stderr.write(bestmodel+" does not have a pdb file!")
        break

    pdb_parser=PDBParser()
    if pdb_path.endswith('.gz'):
        structure=pdb_parser.get_structure(bestmodel, GzipFile(pdb_path))
    else:
        structure=pdb_parser.get_structure(bestmodel, file(pdb_path))

    model = structure[0]
    dssp = DSSP(model, pdb_path, dssp="dssp", acc_array="Wilke")

    buried = []
    exposed = []

    for t in dssp.property_list:
        aa  = t[1]   # one letter codes
        ss  = t[2]   # secondary structure
        rsa = t[3]   # relative surface accessibility

        if rsa >= RSA_EXPOSED_THRESH:
            exposed.append(aa)
        else:
            buried.append(aa)

    exposed = Counter(exposed)
    buried = Counter(buried)

    exposed_total = sum(exposed.values())
    
    exposed_pos_charged = sum([exposed[i] for i in 'KRH'])
    fraction_exposed_pos_charged = float(exposed_pos_charged) / exposed_total

    exposed_neg_charged = sum([exposed[i] for i in 'DE'])
    fraction_exposed_neg_charged = float(exposed_neg_charged) / exposed_total

    exposed_polar = sum([exposed[i] for i in 'QNSTC'])
    fraction_exposed_polar = float(exposed_polar) / exposed_total

    exposed_hydrophobic = sum([exposed[i] for i in 'AVLIMFYWPG'])
    fraction_exposed_hydrophobic = float(exposed_hydrophobic) / exposed_total
    
    query = "UPDATE AA_SURFACE_ACC SET POS_CHARGED='{}', NEG_CHARGED='{}', POLAR='{}', HYDROPHOBIC='{}' WHERE BEST_MODEL='{}'".format(\
    fraction_exposed_pos_charged, fraction_exposed_neg_charged, fraction_exposed_polar, fraction_exposed_hydrophobic, bestmodel)
    c.execute(query)
    
conn.commit()
conn.close()
    

