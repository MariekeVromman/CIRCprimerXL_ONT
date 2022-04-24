#!/usr/bin/python3


import argparse
import os

parser = argparse.ArgumentParser(description='give arguments to filter script')
parser.add_argument("-A", nargs = 1, required=True, help='input file with circRNA + ID')
parser.add_argument("-P", nargs = 1, required=True, help='intput file with all primers')
parser.add_argument('-l', nargs=1, required=True, help='the nr of nucleotides surrounding the BSJ at each side')

args = parser.parse_args()

all_circ = open(args.A[0])
length = int(args.l[0])

# get general info circRNA

circRNA = all_circ.read()

chrom = circRNA.split()[0]
start = int(circRNA.split()[1])
end = int(circRNA.split()[2])
circ_ID = circRNA.split()[3]


# dict of primer annotaiton

anno = open('out_anno.txt', 'r')
an_dict = {}
an_dict_type = {}
for line in anno:
	an_type = line.split('\t')[9]
	exon = line.split('\t')[8]
	primer_id = line.split('\t')[3]
	an_dict[primer_id] = exon
	an_dict_type[primer_id] = an_type

# filter all primers and write to file

all_primers = open(args.P[0])
primer_file = open("all_primers/filtered_primers_" + circ_ID + '_' + chrom + '_' + str(start) + '_' + str(end) + "_.txt", "a")


if os.path.getsize(args.P[0]) == 0:
	design = 0

primer_n = 0

for primer in all_primers:

	# if all tests succeded => primer pair passed filters
	filter_str = "PASS_"


	# add annotation primers
	fw_type = an_dict_type[str(primer_n)+'_F']
	fw_exon = an_dict[str(primer_n)+'_F']

	rv_type = an_dict_type[str(primer_n)+'_R']
	rv_exon = an_dict[str(primer_n)+'_R']


	if fw_type == '.':
		fw_type = 'no_match'
	else:
		if int(fw_type) > 0:
			fw_type = 'intronic/intergenic'
		else:
			fw_type = 'exonic'

	if rv_type == ".":
		rev_type = 'no_match'
	else:
		if int(rv_type) > 0:
			rv_type = 'intronic/intergenic'
		else:
			rv_type = 'exonic'

	exons = 'diff'
	if fw_exon == rv_exon:
		exons = 'same'

	primer_file.write(primer.rstrip() + "\t" + fw_type + "\t" + fw_exon + "\t" + rv_type + "\t" + rv_exon + "\t" + exons + '\t' + filter_str[0:len(filter_str)-1] + "\n")

	primer_n += 1

primer_file.close()
all_primers.close()

