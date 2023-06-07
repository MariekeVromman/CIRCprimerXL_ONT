#!/usr/bin/python3


import argparse
import os

parser = argparse.ArgumentParser(description='give arguments to filter script')
parser.add_argument("-A", nargs = 1, required=True, help='input file with circRNA + ID')
parser.add_argument("-P", nargs = 1, required=True, help='intput file with all primers')

args = parser.parse_args()

# get general info circRNA

circRNA = open(args.A[0]).readline()

chrom = circRNA.split()[0]
start = int(circRNA.split()[1])
end = int(circRNA.split()[2])
circ_ID = circRNA.split()[3]
strand = circRNA.split()[4]

# dict of primer annotation

anno = open('out_anno.txt', 'r')
an_dict = {}
an_dict_type = {}
for line in anno:
	an_type = line.split('\t')[11] # this variable is 0 when the features overlap, or the distance between the features when they don't overlap
	exon = line.split('\t')[8]
	primer_id = line.split('\t')[3]
	an_dict[primer_id] = exon
	an_dict_type[primer_id] = an_type

# get BSJ annotation
bsj_type = an_dict_type['BSJ']

if bsj_type == '.':
	bsj_type = 'no_match'
else:
	if int(bsj_type) > 0:
		bsj_type = 'intronic/intergenic'
	else:
		bsj_type = 'exonic' + '_' + an_dict['BSJ']

# filter all primers and write to file
all_primers = open(args.P[0])
primer_file = open("all_primers/filtered_primers_" + circ_ID + '_' + chrom + '_' + str(start) + '_' + str(end) + "_.txt", "a") # file containing all primers
selected_primer = open('selected_primers_' + circ_ID + '.txt', "w") # file containing selected primer pair
primer_found = "no"


if os.path.getsize(args.P[0]) == 0:
	selected_primer.write(circ_ID + '\t' + chrom + '\t' + str(start) + '\t' + str(end) + '\t' + strand + "\tprimer3 was not able to design primers for this circRNA; try less strict settings\n")
	primer_file.write(circ_ID + '\t' + chrom + '\t' + str(start) + '\t' + str(end) + '\t' + strand  + "\tprimer3 was not able to design primers for this circRNA; try less strict settings\n")
	primer_found = 'yes' # this needs to be changed even if it's not true so that circ is not reported twice

else: 
	primer_n = 0

	for primer in all_primers:

		# add annotation primers
		fw_type = an_dict_type[str(primer_n)+'_F']
		rv_type = an_dict_type[str(primer_n)+'_R']

		if fw_type == '.':
			fw_type = 'no_match'
		else:
			if int(fw_type) > 0: # this variable is 0 when the features overlap, or the distance between the features when they don't overlap
				fw_type = 'intronic/intergenic'
			else:
				fw_type = 'exonic' + '_' + an_dict[str(primer_n)+'_F']

		if rv_type == ".":
			rv_type = 'no_match'
		else:
			if int(rv_type) > 0: # this variable is 0 when the features overlap, or the distance between the features when they don't overlap
				rv_type = 'intronic/intergenic'
			else:
				rv_type = 'exonic' + '_' + an_dict[str(primer_n)+'_R']


		# add column to see if annotations are the same
		exons = 'diff'
		# if they are the same
		if fw_type == rv_type:
			exons = 'same'
		# but if they are both unknown
		if fw_type == 'no_match' or rv_type == "no_match":
			exons = 'no_match'

		# also compare to BSJ annotations

		bsj_same = 'NA'
		if (exons == 'same') & (fw_type == bsj_type):
			bsj_same = 'same'

		primer_file.write(primer.rstrip() + "\t" + fw_type + "\t" + rv_type + "\t" + bsj_type + "\t" + exons + "\t" + bsj_same + "\n")


		# filter best primer pair and save in new output file
		if (primer_found == "no") & (exons == "same"):
			selected_primer.write(primer.rstrip() + "\t" + fw_type + "\t" + rv_type + "\t" + bsj_type + "\t" + exons + "\t" + bsj_same + "\n")
			primer_found = 'yes'

		primer_n += 1

if primer_found == 'no':
	selected_primer.write(circ_ID + '\t' + chrom + '\t' + str(start) + '\t' + str(end) + '\t' + strand + "\tno primer pair could be found with the FWD and REV primer in the same exon\n")

primer_file.close()
all_primers.close()
selected_primer.close()
