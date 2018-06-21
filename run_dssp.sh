#!/bin/bash
# Directory names should NOT end with a '/'
DB=$1   # Path to db
DSSP_OUTDIR=$2   # Path to dir that will store dssp output
MODPIPE_RESULT_DIR=$3   # Path to modpipe results dir
DSSP_PROGRAM=$4   # Path to dssp program
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

