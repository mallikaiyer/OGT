#!/bin/bash
# Diretory names should NOT end with a '/'
DB=$1   # /home/iyer/research/sqlite3/OGT.test.db
DSSP_OUTDIR=$2   # /home/iyer/research/modtest/outputs/dssp_out_rerun 
MODPIPE_RESULT_DIR=$3   # /home/iyer/research/modtest/modpipe_result
DSSP_PROGRAM=$4   # /home/software/bin/dssp
cd ${DSSP_OUTDIR}
for i in `sqlite3 ${DB} "SELECT UNQ_ID, BEST_MODEL FROM BEST_MODELS;"`
do
	UNQ_ID=${i:0:40}
	BEST_MODEL=${i:41:73}
	echo "Working on UNQ_ID: ${UNQ_ID}, BEST_MODEL: ${BEST_MODEL}"
	if [ -f ${BEST_MODEL}.dssp ]; then
		echo "DSSP file for ${BEST_MODEL} already exists."
	else
		${DSSP_PROGRAM} -i ${MODPIPE_RESULT_DIR}/${UNQ_ID:0:3}/${UNQ_ID}/models/${BEST_MODEL}.pdb* -o ${BEST_MODEL}.dssp
		exitcode=$?
		echo Done working on ${BEST_MODEL} with status ${exitcode} @ `date`
	fi
done

