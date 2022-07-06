#!/usr/bin/python3

# # import all libraries
import os
from Bio import Entrez, SeqIO
import argparse

# # get info on BSJ seq
parser = argparse.ArgumentParser(description='give arguments to main primer_xc script')
parser.add_argument('-n', nargs=1, required=True, help='the nr of nucleotides surrounding the BSJ at each side', metavar='length')
parser.add_argument('-i', nargs=1, required=True, help='input circRNA bed file, 0-based')
parser.add_argument('-p', nargs=1, required=True, help='nr of primers')
parser.add_argument('-z', nargs=1, required=True, help='nr of difference between primers')
parser.add_argument('-m', nargs=1, required=True, help='mail address')

parser.add_argument('-a', nargs=1, required=True, help='min TM')
parser.add_argument('-b', nargs=1, required=True, help='max TM')
parser.add_argument('-c', nargs=1, required=True, help='opt TM')
parser.add_argument('-d', nargs=1, required=True, help='TM diff')
parser.add_argument('-e', nargs=1, required=True, help='min GC')
parser.add_argument('-f', nargs=1, required=True, help='max GC')
parser.add_argument('-g', nargs=1, required=True, help='opt GC')
parser.add_argument('-j', nargs=1, required=True, help='min amp length')
parser.add_argument('-k', nargs=1, required=True, help='max amp length')

args = parser.parse_args()
length = int(args.n[0])
input_bed = open(args.i[0])
nr = args.p[0]
diff = args.z[0]
mail = args.m[0]


# # retrieve one side of BSJ (depending on strand info)

# ## change chrom nr to GI from NCBI Revesion History GRCh38.p13 Primary Assembly
chr_to_GI = {"chr1":"568815597","chr2":"568815596", "chr3":"568815595", "chr4": "568815594", "chr5":"568815593", 'chr6':"568815592", "chr7":"568815591", "chr8":'568815590', "chr9":"568815589", "chr10":"568815588", "chr11":"568815587", "chr12":"568815586", "chr13":"568815585", "chr14":"568815584", "chr15":"568815583", "chr16":"568815582", "chr17":"568815581", "chr18":"568815580", "chr19":"568815579", 'chr20':"568815578", "chr21":"568815577", "chr22":"568815576", "chrX":"568815575", "chrY":"568815574"}

# # read first (and only) line of input bed file
circRNA = input_bed.read()

# retrieve chrom start end info
chrom = circRNA.split()[0]
start = int(circRNA.split()[1]) + 1 # change to 1-based system
end = int(circRNA.split()[2])
circ_ID = circRNA.split()[3]
circ_strand = circRNA.split()[4]

# change chr nr to GI nr
chrom_GI = chr_to_GI[chrom]

# ## retrieve right side of BSJ

from Bio import Entrez, SeqIO
Entrez.email = "marieke.vromman@ugent.be"

if circ_strand == "+":
    handle = Entrez.efetch(db="nucleotide", 
                       id=chrom_GI, 
                       rettype="fasta", 
                       strand=1, 
                       seq_start=start+30, 
                       seq_stop=start+length-1)
elif circ_strand == '-':
    handle = Entrez.efetch(db="nucleotide", 
                       id=chrom_GI, 
                       rettype="fasta", 
                       strand=2, 
                       seq_start=end-30, 
                       seq_stop=end-length)
else:
    raise SystemExit('{0} is not a valid strand for circ {1}'.format(circ_strand, circRNA))


record_right = SeqIO.read(handle, "fasta")
handle.close()


# ## change into string
sequence = str(record_right.seq)

# ## make a txt file for primer design (sequence based)

output = open("input_primer3_" + circ_ID + ".txt", "w")
output.write("SEQUENCE_ID=" + circ_ID + "_" + chrom + "_" + str(start-1) + "_" + str(end) + '_' + circ_strand + "\n")
output.write("SEQUENCE_TEMPLATE=" + sequence + "\n")
#output.write("SEQUENCE_TARGET=" + str(30) + ",1\nPRIMER_NUM_RETURN=" + nr + "\n")
output.write("PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE=" + diff + '\n')
output.write("PRIMER_PRODUCT_SIZE_RANGE=" + args.j[0] + '-' + args.k[0] + '\n')
output.write("PRIMER_MIN_TM=" + args.a[0] + '\n')
output.write("PRIMER_MAX_TM=" + args.b[0]+ '\n')
output.write("PRIMER_OPT_TM=" + args.c[0]+ '\n')
output.write("PRIMER_PAIR_MAX_DIFF_TM=" + args.d[0]+ '\n')
output.write("PRIMER_MIN_GC="  + args.e[0]+ '\n')
output.write("PRIMER_MAX_GC=" + args.f[0]+ '\n')
output.write("PRIMER_OPT_GC_PERCENT="  + args.g[0]+ '\n')
output.write("=\n")
output.close()


# ## make file containing general circ info

output = open("circ_info_" + circ_ID + ".bed", 'w')
output.write(circRNA)
output.close()


