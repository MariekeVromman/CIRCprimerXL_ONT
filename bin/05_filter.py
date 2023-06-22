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
an_dict_all = {}
for line in anno:
	an_type = line.split('\t')[11] # this variable is 0 when the features overlap, or the distance between the features when they don't overlap
	exon = line.split('\t')[8]
	primer_id = line.split('\t')[3]
	an_dict[primer_id] = exon   # {'1_F' : 'ENSG00000166532,RIMKLB,ENST00000619374,1', 'BSJ' : 'ENSG00000166532,RIMKLB,ENST00000619374,1'}
	an_dict_type[primer_id] = an_type # {'1_F' : '0', 'BSJ' : '0'}
	an_dict_all[primer_id] = line # dict that has all info, we might need later

print(an_dict)
print(an_dict_type)

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

	for primer in all_primers:

		primer_n = primer.split('\t')[5]

		# add annotation primers
		fw_type = an_dict_type[str(primer_n)+'_F']
		rv_type = an_dict_type[str(primer_n)+'_R']

		# get information about overlap to check if overlap is complete or not
		# chr12	8713810	8713811	BSJ	+	chr12	8713742	8714041	ENSG00000166532,RIMKLB,ENST00000619374,1	0	+	0

		fw_primer_s = an_dict_all[str(primer_n)+'_F'].split('\t')[1]
		fw_primer_e = an_dict_all[str(primer_n)+'_F'].split('\t')[2]
		fw_match_s = an_dict_all[str(primer_n)+'_F'].split('\t')[6]
		fw_match_e = an_dict_all[str(primer_n)+'_F'].split('\t')[7]

		overlap_lab_fw = ''

		if (fw_primer_s < fw_match_s) or (fw_primer_e > fw_match_e):
			overlap_lab_fw = 'partially_'

		rv_primer_s = an_dict_all[str(primer_n)+'_R'].split('\t')[1]
		rv_primer_e = an_dict_all[str(primer_n)+'_R'].split('\t')[2]
		rv_match_s = an_dict_all[str(primer_n)+'_R'].split('\t')[6]
		rv_match_e = an_dict_all[str(primer_n)+'_R'].split('\t')[7]

		overlap_lab_rv = ''

		if (rv_primer_s < rv_match_s) or (rv_primer_e > rv_match_e):
			overlap_lab_rv = 'partially_'

		# determine match
		if fw_type == '.':
			fw_type = 'no_match'
		else:
			if int(fw_type) > 0: # this variable is 0 when the features overlap, or the distance between the features when they don't overlap
				fw_type = 'intronic/intergenic'
			else:
				fw_type = overlap_lab_fw + 'exonic' + '_' + an_dict[str(primer_n)+'_F']

		if rv_type == ".":
			rv_type = 'no_match'
		else:
			if int(rv_type) > 0: # this variable is 0 when the features overlap, or the distance between the features when they don't overlap
				rv_type = 'intronic/intergenic'
			else:
				rv_type = overlap_lab_rv + 'exonic' + '_' + an_dict[str(primer_n)+'_R']


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
		if (primer_found == "no") & (exons == "same") & (overlap_lab_fw == "") & (overlap_lab_rv == ""):
			selected_primer.write(primer.rstrip() + "\t" + fw_type + "\t" + rv_type + "\t" + bsj_type + "\t" + exons + "\t" + bsj_same + "\n")
			primer_found = 'yes'

if primer_found == 'no':
	selected_primer.write(circ_ID + '\t' + chrom + '\t' + str(start) + '\t' + str(end) + '\t' + strand + "\tno primer pair could be found with the FWD and REV primer in the same exon\n")

primer_file.close()
all_primers.close()
selected_primer.close()
