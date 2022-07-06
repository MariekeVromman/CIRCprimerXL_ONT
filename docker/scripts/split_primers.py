#!/usr/bin/python3

rev_comp_lib = {'A': 'T', 'T':'A', 'G':'C', 'C':'G'}

def rev_comp(seq):
	rev_seq = ''
	for nt in seq:
		rev_seq = rev_seq + rev_comp_lib[nt]
	rev_seq = rev_seq[::-1]
	return rev_seq

import argparse

parser = argparse.ArgumentParser(description='give arguments to main primer_xc script')
parser.add_argument('-i', nargs=1, required=True, help='input primer file')
args = parser.parse_args()
input_primers = args.i[0]

primer_in = open(input_primers)

primers = {}

circ_info_keys = ("SEQUENCE_ID", "SEQUENCE_TEMPLATE", "SEQUENCE_TARGET")

# read all info into dictionary
for line in primer_in:
	key, value = line.split("=")
	value = value.rstrip()
	primers[key] = value

template = primers["SEQUENCE_TEMPLATE"]
circRNA = primers["SEQUENCE_ID"]
circ_ID, chrom, start, end, strand = circRNA.split("_")
nr_p_out = primers["PRIMER_LEFT_NUM_RETURNED"]

primer_in.close()

# extra loop to copy everything to a dir containing primer3 details
primer_in = open(input_primers)
primer_extra = open("output_primer3_" + circ_ID + '.txt', 'a')
for line in primer_in:
	primer_extra.write(line)
primer_in.close()
primer_extra.close()


# make general file with list primers (cleans up primer3 output)
all_primers = open("all_primers_" + circ_ID + ".txt", 'w')
all_primers_dict = {}

# and simlyltaniously make bed file that can be used to annotate the circRNA
bed_in = open('bed_in.txt', 'w')


for primer_index in range(int(nr_p_out)):

	FWD = primers[("PRIMER_LEFT_" + str(primer_index) + "_SEQUENCE")]
	REV = primers[("PRIMER_RIGHT_" + str(primer_index) + "_SEQUENCE")]

	FWD_pos, FWD_len = primers['PRIMER_LEFT_'+ str(primer_index)].split(",")
	REV_pos, REV_len = primers['PRIMER_RIGHT_'+ str(primer_index)].split(",")
	amplicon = template[int(FWD_pos):int(REV_pos) + 1]

	PRIMER_LEFT_TM = primers[("PRIMER_LEFT_" + str(primer_index) + "_TM")]
	PRIMER_RIGHT_TM = primers[("PRIMER_RIGHT_" + str(primer_index) + "_TM")]
	PRIMER_LEFT_GC_PERCENT = primers[("PRIMER_LEFT_" + str(primer_index) + "_GC_PERCENT")]
	PRIMER_RIGHT_GC_PERCENT = primers[("PRIMER_RIGHT_" + str(primer_index) + "_GC_PERCENT")]

	FWD_RC = rev_comp(FWD)
	REV_RC = rev_comp(REV)

	## make file for bedtools anno

	# FWD primer
	bed_in.write(chrom + "\t" + str(int(start)+30+int(FWD_pos)) + "\t" + str(int(start)+30+int(FWD_pos)+int(FWD_len)) + '\t'+ str(primer_index) + '_F' + '\t' + strand +'\n')
	# REV primer
	bed_in.write(chrom + "\t" + str(int(start)+30+int(REV_pos)-int(REV_len)) + "\t" + str(int(start)+30+int(REV_pos)) + '\t'+ str(primer_index) + '_R' + '\t' + strand + '\n')


	# general primer file (for filtering), first put in dict, will be sorted (see below)
	all_primers_dict[circ_ID + "\t" + chrom + "\t" + start + "\t" + end + '\t' + strand + '\t' + str(primer_index) + '\t' + FWD + '\t' + FWD_RC + '\t' + REV + '\t' + 
	REV_RC + '\t' + FWD_pos + '\t' + FWD_len + '\t' + REV_pos +'\t' + REV_len + '\t' + PRIMER_LEFT_TM + '\t' + PRIMER_RIGHT_TM + '\t' + 
	PRIMER_LEFT_GC_PERCENT + '\t' + PRIMER_RIGHT_GC_PERCENT + '\t' + amplicon + '\t' + '\n'] = len(amplicon)


# also add BSJ to bed file
bed_in.write(chrom + "\t" + str(start) + "\t" + str(int(start)+1) + '\t'+ 'BSJ' + '\t' + strand +'\n')

# sort primers according to amp size (smallest is best) and then print to all_primers
all_primers_sorted = {k: v for k, v in sorted(all_primers_dict.items(), key=lambda item: item[1])}

for primer in all_primers_sorted:
	all_primers.write(primer)


all_primers.close()
bed_in.close()



