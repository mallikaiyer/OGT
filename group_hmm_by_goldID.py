#!/usr/bin/env python2
"""
IN: csv list of gold ID to Genbank ID mappings
    hmmsearch output files run on each Genbank ID 
    
OUT: produce a single hmmsearch output file with all hits for all HMMs
     above TC per GoldID
     produce a single fasta file with <= 10 records, with BEST hit above TC
     kept for each hmm
"""

import csv
import os
import sys
import gzip
from Bio import SeqIO

#

GOLD_IDS_TO_GENBANK = "/home/chris/OGT/modpipe/goldID_genbankID.csv"

## e.g. /home/iyer/research/modtest/outputs/hmmcheck/by_genbank/GENBANK/HMM.out
hmm_output_dir = "/home/iyer/research/modtest/outputs/hmmcheck"
## eg.  /home/iyer/research/genomes_ncbi/genomes/GENBANK/*_protein.faa*
genbank_dirs   = "/home/iyer/research/genomes_ncbi/genomes"
fasta_output_dir = "/home/iyer/research/modtest/proteomes/proteomes_of_ten"

# first keep a dictionary of lists to map goldID to genbankID
gold2genbank = {}
for row in csv.reader(file(GOLD_IDS_TO_GENBANK)):
    goldID, genbankID = row
    if goldID not in gold2genbank:
        gold2genbank[goldID] = set()
    gold2genbank[goldID].add(genbankID)

for goldID, genbankset in gold2genbank.iteritems():
    print goldID
    # concatenate all HMM.out files from all genbank entries into one HMM file.
    all_genbank_hmm_hits_fn = os.path.join(hmm_output_dir, "by_goldid", goldID+"_hmmhits.out") 
    if os.path.isfile(all_genbank_hmm_hits_fn):
        os.remove(all_genbank_hmm_hits_fn)  # erase file if exists
    for genbankID in genbankset:
        print " "+genbankID
        os.system("""egrep -v '^#' {0} | awk '{{if($1!="#")print "{1} "$0}}' >> {2}""".format(os.path.join(hmm_output_dir, "by_genbank", genbankID, "*HMM.out"), genbankID, all_genbank_hmm_hits_fn )) 
    
    # go through concatenated file, get best hit for each hmm across all genbank files
    hmm_hits = {}
    with open(all_genbank_hmm_hits_fn, 'r') as all_genbank_hmm_hits:
        for line in all_genbank_hmm_hits:
            if line.startswith('#'):
                continue
            #print line
            items = line.split(None, 19)  # hmm output is space separated into 19 fields, and we have added one more at beginning (genbankID), with last being description
            #print items
            hmm = items[4]
            evalue = items[5]
            genbankID = items[0]
            accession = items[1]  # the actual protein sequence hit
            print "     "+accession
            hmm_hits.setdefault(hmm, [])
            hmm_hits[hmm].append((evalue, genbankID, accession))
    # now have a dict with tuples keyed on hmm.  Go through each hmm and pick best with sorting.
    # write out the BEST hmm hit for each of ten hmms to a new fasta file
    with open(os.path.join(fasta_output_dir, goldID+"_proteome_ten.faa"), 'w') as output_fasta:
        for hmm in sorted(hmm_hits):
            # sort by evalue, take lowest evalue (pos 0) and we have best hit
            evalue, genbankID, accession = sorted(hmm_hits[hmm])[0]
            # go get that specific fasta sequence, and write to output file with modified header
            faa_path = None
            for candidate_filename in os.listdir(os.path.join(genbank_dirs, genbankID)):
                if candidate_filename.endswith(("_protein.faa.gz", "_protein.faa")):
                    faa_path = os.path.join(genbank_dirs, genbankID, candidate_filename)
                    break
                else:
                    continue
            if faa_path is None:
                raise ValueError, "No proteome file for genbank ID {}".format(genbankID)
            elif faa_path.endswith(".gz"):
                faa_file_handle = gzip.open(faa_path, 'r')
            else:
                faa_file_handle = open(faa_path, 'r')
            # dumb, but scan through file for that id
            #print "     "+accession 
            for record in SeqIO.parse(faa_file_handle, "fasta"):
                #print record
                if record.id == accession:
                    print record.description
                    record.description = "{} hmmhit={} hmmevalue={} {}".format(record.id, hmm, evalue, record.description.split(None, 1)[1])
                    output_fasta.write(record.format("fasta"))
            faa_file_handle.close()

    break #DEBUG
