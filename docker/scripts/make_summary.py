#!/usr/bin/python3

import os

# get list of circ
all_circ = open('filtered_primers.txt')
all_circ_ls = []

primer_found = 0
primer3_fail = 0
exon_fail = 0

# skip the first line (header)
all_circ.readline()

# add circ level info
for circ in all_circ:
    all_circ_ls.append(circ.split()[0])

    if len(circ.split('\t')) > 7:
        primer_found += 1
    elif circ.split()[5] == 'primer3':
        primer3_fail += 1
    elif circ.split()[5] == 'no':
        exon_fail +=1

nr_circ = len(all_circ_ls)

summary = open('run_summary.txt', 'w')
summary.write('in this run, primers were designed for {0} circRNAs\n\tfor {1} circRNAs primer could be designed\n\tfor {2} circRNAs no primers could be designed with the same FWD and REV exon annotation\n\tfor {3} circRNAs no primers could be designed by primer3'.format(nr_circ, primer_found, exon_fail, primer3_fail))
summary.close()
