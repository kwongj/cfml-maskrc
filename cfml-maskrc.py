#!/usr/bin/env python
# Script by Jason Kwong
# Script to mask recombination from CFML output

# Usage
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(
	formatter_class=RawTextHelpFormatter,
	description='Script to mask recombination from ClonalFrameML output',
	usage='\n  %(prog)s --aln FASTA --out OUTFILE --symbol ? <CFML_PREFIX>')
parser.add_argument('prefix', metavar='CFML_PREFIX', help='prefix used for CFML output files')
parser.add_argument('--aln', metavar='FASTA', required=True, help='multiFASTA alignment used as input for CFML')
parser.add_argument('--out', metavar='OUTFILE', default='maskrc.aln', help='output file for masked alignment (default="maskrc.aln")')
parser.add_argument('--symbol', metavar='?', default='?', help='symbol to use for masking (default="?")')
parser.add_argument('--version', action='version', version=
	'=====================================\n'
	'%(prog)s v0.1\n'
	'Updated 2-Oct-2015 by Jason Kwong\n'
	'Dependencies: Python 2.x, BioPython, ete2\n'
	'=====================================')
args = parser.parse_args()
cfmlPREFIX = str(args.prefix)
cfmlTREE = cfmlPREFIX + ".labelled_tree.newick"
cfmlRECOMB = cfmlPREFIX + ".importation_status.txt"
outfile = args.out
symbol = args.symbol

import csv
from collections import defaultdict
from ete2 import Tree
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord

# Import CFML files and identify recombinant regions
d = defaultdict(list)
seqLIST = []
t = Tree(cfmlTREE, format=1)
with open(cfmlRECOMB) as csvfile:
	RCseqs = csv.reader(csvfile, delimiter='\t')
	next(csvfile)
	for row in RCseqs:
		seq = row[0]
		RCstart = row[1]
		RCstop = row[2]
		node = t.search_nodes(name=seq)[0]
		if node.is_leaf():
			leaf_names = [node.name]
		else:
			leaves = node.get_leaves()
			leaf_names = [leaf.name for leaf in leaves]
		for l in leaf_names:
			seqLIST.append([l, [RCstart, RCstop]])

# Build dictionary of recombinant regions
for k,v in seqLIST:
	d[k].append(v)

# Mask recombination in sequences
seqALN = []
for record in SeqIO.parse(args.aln, 'fasta'):
	print 'Reading %s ... ' % record.id
	regions = d.get(record.id, None)
	newrec = MutableSeq(str(record.seq))
	if regions:
		for a in regions:
			start = int(a[0]) - 1
			end = int(a[1])		
			lenMASK = end - start
			newrec[start:end] = (symbol)*lenMASK
	seqALN.append(SeqRecord(Seq(str(newrec)), record.id, description=''))

# Write masked alignment to file
print 'Writing masked alignment to %s ... ' % outfile
SeqIO.write(seqALN, outfile, 'fasta')
