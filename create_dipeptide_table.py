#!/usr/bin/env python2

'''
Creates table to store all dipeptide frequencies. This script also populates the table with goldIDs and OGTs.
'''
import sys
import sqlite3

DB = sys.argv[1]
conn = sqlite3.connect(DB)
c = conn.cursor()

query = "CREATE TABLE DIPEPTIDE_FREQUENCIES (GOLD_ID VARCHAR(15) PRIMARY KEY, OGT REAL" 

aa_list = ['G', 'A', 'I', 'L', 'V', 'M', 'C', 'D', 'E', 'N', 'Q', 'K', 'R', 'H', 'Y', 'W', 'F', 'S', 'T', 'P']

for i in aa_list:
	for j in aa_list:
		dipeptide = i + j
		query = query + ',' + ' ' + 'F_{}'.format(dipeptide) + ' ' + 'REAL'

query = query + ');'
c.execute(query)

query = "INSERT INTO DIPEPTIDE_FREQUENCIES (GOLD_ID, OGT) SELECT DISTINCT(GOLD_ID), TEMPERATURE_OPTIMUM_CLEAN FROM temp_acc WHERE TEMPERATURE_OPTIMUM_CLEAN!='';"
c.execute(query)
conn.commit()
conn.close()
