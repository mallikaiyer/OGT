#!/usr/bin/env python2

'''
This script selects UNQ_ID from the GOLD_PROTEINS table and finds the model 
with the lowest DOPE score. This best model is added to the BEST_MODELS table.
Input: (1) Path to dir containing modpipe result.
(2) Path to db.
'''

import os
from os import path
import gzip
import sqlite3
import sys


modpipe_result_dir = sys.argv[1] # Path to dir containing modpipe results. 
db = sys.argv[2] # Path to db
conn = sqlite3.connect(db) 
c = conn.cursor()

query = "SELECT DISTINCT(UNQ_ID) FROM GOLD_PROTEINS;"
c.execute(query)
unqIDs = [i[0] for i in c.fetchall()]
for unqID in unqIDs:
    print "Working on {}".format(unqID)
    score_dict = {}
    models_dir = os.path.join(modpipe_result_dir, unqID[0:3], unqID, 'models')
    if os.path.isdir(models_dir):
        for file in os.listdir(models_dir): 
            filename = os.path.join(models_dir, file)
            if filename.endswith(".gz"):
                f = gzip.open(filename, 'rb')
            else:
                f = open(filename, 'r')
            
            for j, line in enumerate(f):
                if j == 6: # DOPE score REMARK should be the 7th line in the file.
                    l = line.split()
                    if len(l) > 5:   # Simple check to make sure that l is the line with the DOPE score REMARK.
                        print "This is not the DOPE score remark for {}".format(filename)
                    else:
                        modelID = file[0:32]
                        score = float(l[4]) # The fifth element in the line is the DOPE score.
                        score_dict[modelID] = score

    else:
	   print "{} does not exist.".format(models_dir)

    for key in score_dict:
        if score_dict[key] == min(score_dict.values()):
            conn.execute("INSERT INTO BEST_MODELS (UNQ_ID, BEST_MODEL, DOPE_SCORE) VALUES ('{}', '{}', '{}');".format(unqID, key, score_dict[key]))    

conn.commit()
conn.close()
            
print "Done."
