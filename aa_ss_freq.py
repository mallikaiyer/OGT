#!/usr/bin/env python2
'''
Input: 
(1) Path to db
(2) Path to dir with DSSP output files.

This script pulls the BEST_MODELs from the BEST_MODELS table. For each best model, the dssp output file
is parsed into a dict by the PDB.make_dssp_dict function. This dict is then used to calculate the frac
of positively-charged, negatively-charged, polar and hydrophobic amino acids in helices, beta sheets and 
loops. The AA_HELIES, AA_BETA and AA_LOOPS tables are updated with these fractions.

Run as follows: ./aa_ss_freq_publish.py /path/to/db/ /path/to/dssp/output/dir/

'''

import sys
import sqlite3
import os
from Bio import PDB
from collections import Counter



if __name__ == '__main__':


    db = sys.argv[1] # Path to db.
    dssp_dir = sys.argv[2] # Path to dir with DSSP output files.
    conn = sqlite3.connect(db) 
    c = conn.cursor()

    query = "SELECT BEST_MODEL FROM BEST_MODELS;"
    c.execute(query)
    bestModels = [str(i[0]) for i in c.fetchall()]
    for bestModel in bestModels:
        dssp_file = os.path.join(dssp_dir, bestModel+'.dssp')
        d, keys = PDB.make_dssp_dict(dssp_file)
        ss_list    = []
        aa_list    = []
        cys_set = {"C", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"}
        helix_set = {"H", "G", "I"}
        beta_set = {"B", "E"} 
        loop_set = {"T", "S"} 

        helix_res  = []
        beta_res = []
        loop_res = []

        for k in keys:
            chain, res = k
            res_number = res[1]
            dssp_tuple = d[k]
            aa = dssp_tuple[0]
            ss = dssp_tuple[1]
            solvent_acc = dssp_tuple[2]

            # populate lists
            ss_list.append(ss)
            aa_list.append(aa)
            if ss in helix_set:
                helix_res.append(aa)
            elif ss in beta_set:
                beta_res.append(aa)
            elif ss in loop_set:
                loop_res.append(aa)
            elif ss == "-":
                continue
            else:
                print "{} is not a secondary_struc code. (res_number:{}, dssp_file:{})".format(ss, res_number, dssp_file)
            
        ch = Counter(helix_res)
        cb = Counter(beta_res)
        cl = Counter(loop_res)

        if sum(ch.values()) != 0:
            positive_helix = float(ch['K'] + ch['H'] + ch['R']) / sum(ch.values())
            negative_helix = float(ch['D'] + ch['E']) / sum(ch.values())
            Ch = float(sum(ch[x] for x in cys_set)) / sum(ch.values())
            polar_helix = float(ch['Q'] + ch['N'] + ch['S'] + ch['T'] + Ch) / sum(ch.values())
            hydrophobic_helix = float(ch['A'] + ch['F'] + ch['I'] + ch['L'] + ch['M'] + ch['V'] + ch['W'] + ch['Y']) / sum(ch.values())
            p_helix = float(ch['P']) / sum(ch.values())
            g_helix = float(ch['G']) / sum(ch.values())

        else:
            print "No helix residues in best model: {}".format(bestModel)
            positive_helix = 0
            negative_helix = 0
            polar_helix = 0
            hydrophobic_helix = 0
            p_helix = 0
            g_helix = 0

        if sum(cb.values()) != 0:
            positive_beta = float(cb['K'] + cb['H'] + cb['R']) / sum(cb.values())
            negative_beta = float(cb['D'] + cb['E']) / sum(cb.values())
            Cb = float(sum(cb[x] for x in cys_set)) / sum(cb.values())
            polar_beta = float(cb['Q'] + cb['N'] + cb['S'] + cb['T'] + Cb) / sum(cb.values())
            hydrophobic_beta = float(cb['A'] + cb['F'] + cb['I'] + cb['L'] + cb['M'] + cb['V'] + cb['W'] + cb['Y']) / sum(cb.values())
            p_beta = float(cb['P']) / sum(cb.values())
            g_beta = float(cb['G']) / sum(cb.values())

        else:
            print "No beta sheet/strand residues in best model: {}".format(bestModel)
            positive_beta = 0
            negative_beta = 0
            polar_beta = 0
            hydrophobic_beta = 0
            p_beta = 0
            g_beta = 0

        if sum(cl.values()) != 0:
            positive_loops = float(cl['K'] + cl['H'] + cl['R']) / sum(cl.values())
            negative_loops = float(cl['D'] + cl['E']) / sum(cl.values())
            Cl = float(sum(cl[x] for x in cys_set)) / sum(cl.values())
            polar_loops = float(cl['Q'] + cl['N'] + cl['S'] + cl['T'] + Cl) / sum(cl.values())
            hydrophobic_loops = float(cl['A'] + cl['F'] + cl['I'] + cl['L'] + cl['M'] + cl['V'] + cl['W'] + cl['Y']) / sum(cl.values())
            p_loops = float(cl['P']) / sum(cl.values())
            g_loops = float(cl['G']) / sum(cl.values())

        else:
            print "No loop residues in best model: {}".format(bestModel)
            positive_loops = 0
            negative_loops = 0
            polar_loops = 0
            hydrophobic_loops = 0
            g_loops = 0
            p_loops = 0
        
        query = "UPDATE AA_HELICES SET POS_CHARGED={}, NEG_CHARGED={}, POLAR={}, HYDROPHOBIC={}, P={}, G={} WHERE BEST_MODEL='{}';".format(positive_helix, negative_helix, polar_helix, hydrophobic_helix, p_helix, g_helix, bestModel)
        conn.execute(query)
        query = "UPDATE AA_BETA SET POS_CHARGED={}, NEG_CHARGED={}, POLAR={}, HYDROPHOBIC={}, P={}, G={} WHERE BEST_MODEL='{}';".format(positive_beta, negative_beta, polar_beta, hydrophobic_beta, p_beta, g_beta, bestModel)
        conn.execute(query)
        query = "UPDATE AA_LOOPS SET POS_CHARGED={}, NEG_CHARGED={}, POLAR={}, HYDROPHOBIC={}, P={}, G={} WHERE BEST_MODEL='{}';".format(positive_loops, negative_loops, polar_loops, hydrophobic_loops, p_loops, g_loops, bestModel)
        conn.execute(query)
        conn.commit()

    conn.close()
