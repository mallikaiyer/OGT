#!/usr/bin/env python2
'''
Input: Path to db.

This script calculates the fraction of each amino acid group in helices, beta sheets/strands, loops and 
exposed regions by gold ID. It calculates the average of each feature for all the best models of one UNQ_ID.
It then caculates the average of each feature for all UNQ_IDs of one GOLD_ID and updates the GOLD_FEATURES_TABLE
with that value. This is done for all GOLD_IDs.

Run as follows: ./average_aa_ss_surface_features.py /path/to/db/
'''

import sqlite3
import sys

db = sys.argv[1] # Path to db
conn = sqlite3.connect(db)
c = conn.cursor()

def average_feature(struc, aaGroup):

    query = "SELECT DISTINCT(UNQ_ID) FROM BEST_MODELS;"
    c.execute(query)
    unqIDs = [str(i[0]) for i in c.fetchall()]
    unq_bestmodels = {}
    unq_average_feature = {}
    for unqID in unqIDs:
        print 'Working on unqID: ', unqID
        query = "SELECT BEST_MODEL FROM BEST_MODELS WHERE UNQ_ID = '{}';".format(unqID)
        c.execute(query)
        unq_bestmodels[unqID] = [str(i[0]) for i in c.fetchall()]
        unq_average_feature[unqID] = 0
        unq_best_model_features = []
        for bestModel in unq_bestmodels[unqID]:
            print 'Working on best model: ', bestModel
            query = "SELECT {} from AA_{} WHERE BEST_MODEL = '{}';".format(aaGroup, struc, bestModel)
            c.execute(query)
            featureValue = c.fetchone()[0]
            unq_best_model_features.append(featureValue)
            unq_average_feature_value = sum(unq_best_model_features) / float(len(unq_best_model_features))    
            unq_average_feature[unqID] = unq_average_feature_value

    query = "SELECT DISTINCT GOLD_ID FROM GOLD_PROTEINS;" 
    c.execute(query)
    gold_list = [i[0] for i in c.fetchall()]
    for goldID in gold_list:
        print 'Working on goldID: ', goldID
        # query = "SELECT UNQ_ID FROM GOLD_PROTEINS WHERE GOLD_ID = '{}';".format(goldID)
        query = "SELECT UNQ_ID FROM GOLD_PROTEINS WHERE GOLD_ID = '{}' AND UNQ_ID IN (SELECT DISTINCT(UNQ_ID) \
        FROM BEST_MODELS);".format(goldID) # The GOLD_PROTEINS table has some UNQ_IDs that don't have models.
        c.execute(query)
        gold_unq_list = [str(i[0]) for i in c.fetchall()] 
        if len(gold_unq_list) == 0:
            sys.stderr.write('{} has NO modelled proteins! Continuing to next gold ID'.format(goldID))
            continue
        gold_unq_features = [unq_average_feature[unqID] for unqID in gold_unq_list]
        gold_average_feature = sum(gold_unq_features) / float(len(gold_unq_features))

        query = "UPDATE GOLD_FEATURES SET {}_{} = {} WHERE GOLD_ID = '{}';".format(struc, aaGroup, gold_average_feature, goldID)
        c.execute(query)
        conn.commit()


if __name__ == '__main__':
    strucTypes = ['HELICES', 'BETA', 'LOOPS']
    aaGroups = ['POS_CHARGED', 'NEG_CHARGED', 'POLAR', 'HYDROPHOBIC', 'P', 'G']
    for struc in strucTypes:
        for aaGroup in aaGroups:
            average_feature(struc, aaGroup)

    aaGroups_surfaceAcc = ['POS_CHARGED', 'NEG_CHARGED', 'POLAR', 'HYDROPHOBIC']
    for aaGroup_surfaceAcc in aaGroups_surfaceAcc:
        average_feature('SURFACE_ACC', aaGroup_surfaceAcc)

    conn.close()

