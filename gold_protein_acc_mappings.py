#!/usr/bin/env python2

import sys
import os
import sqlite3

if __name__ == '__main__':

	proteomes_dir = sys.argv[1] #  Path to dir containing proteomes of ten proteins each, for each gold ID. '/home/iyer/research/modtest/proteomes/proteomes_of_ten/'
	db = sys.argv[2] #  Path to db. '/home/iyer/research/sqlite3/OGT.test.db'
	conn = sqlite3.connect(db)
	for file in os.listdir(proteomes_dir):
		gold_id = file[0:-17]
		with open(os.path.join(proteomes_dir, file)) as handle:
			for line in handle:
				if line.startswith('>'):
					atoms = line.split()
					prot_acc = atoms[0][1:]
					# print gold_id+","+prot_acc
					query = "INSERT INTO GOLD_PROTEINS (GOLD_ID, PROTEIN_ACC) VALUES ('{}', '{}');".format(gold_id, prot_acc)
					conn.execute(query)
					conn.commit()
				else:
					continue

	conn.close()
