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


# filter all primers and write to file

all_primers = open(args.P[0])
primer_file = open("all_primers/filtered_primers_" + circ_ID + '_' + chrom + '_' + str(start) + '_' + str(end) + "_.txt", "a")


if os.path.getsize(args.P[0]) == 0:
	design = 0


for primer in all_primers:

	# if all tests succeded => primer pair passed filters
	filter_str = "PASS_"


	primer_file.write(primer.rstrip() + "\t" + filter_str[0:len(filter_str)-1] + "\n")

primer_file.close()
all_primers.close()

