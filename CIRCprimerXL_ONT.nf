#!/usr/bin/env nextflow

// set default parameters

params.primer_settings = "$baseDir/assets/primer3plus_settings.txt"
params.chrom_file = "$baseDir/assets/GRCh38/chrom_sizes_GRCh38.txt"
params.input_bed = "example/path"
params.exons = "$baseDir/assets/GRCh38/Homo_sapiens.GRCh38.103.exons.sorted.bed"
params.index_fasta = "$baseDir/assets/GRCh38/index_fastahack"
params.index_fasta_name = "GRCh38.dna.primary_assembly.fa"
params.index_bowtie = "$baseDir/assets/GRCh38/index_bowtie"
params.index_bowtie_name = "GRCh38_dna"


params.primer3_diff = 1
params.primer3_nr = 20
params.min_tm = 60
params.max_tm = 70
params.opt_tm = 65
params.diff_tm = 2
params.min_gc = 40
params.max_gc = 80
params.opt_gc = 50
params.amp_min = 40
params.amp_max = 1000
params.temp_l = 250
params.design = '5_prime'


// change parameters to files
input_bed = file(params.input_bed)
chrom_file = file(params.chrom_file)
exons = file(params.exons)

// print run info

log.info """\
==============================================
CIRCprimerXL pipeline ONT
==============================================
OncoRNALab - Marieke Vromman
https://github.com/MariekeVromman/CIRCprimerXL_ONT
https://hub.docker.com/repository/docker/oncornalab/circprimerxl_ont
==============================================
your input file : ${params.input_bed}
your output directory : ${params.output_dir}
"""

// define and run each process

process split_circRNAs {

	input:
	path 'input_bed_handle' from input_bed
	path 'chrom_file_handle' from chrom_file

	output:
	path 'circ*' into ind_circ_file

	"""
	01_validate_bed.py -i $input_bed_handle -c $chrom_file_handle
	02_split_circRNAs.py -i $input_bed_handle
	"""

}


process get_seq {
	maxForks 1 // this parameter is necessary as Entrez Bio only allows a limited nr of request per mail address in th get_circ_seq.py script
	input:
	file ind_circ_file_handle from ind_circ_file.flatten()
	val 'diff' from params.primer3_diff
	val 'nr' from params.primer3_nr
	val 'length' from params.temp_l

	output:
	path 'input_primer3*' into in_primer3
	path 'circ_info*' into circ_info

	"""
	03_get_circ_seq.py -n $length -i $ind_circ_file_handle -z $diff -p $nr -m 'marieke.vromman@ugent.be' -a $params.min_tm -b $params.max_tm -c $params.opt_tm -d $params.diff_tm -e $params.min_gc -f $params.max_gc -g $params.opt_gc -j $params.amp_min -k $params.amp_max -l $params.design
	"""
}



process get_primers {

	publishDir "$params.output_dir/primer3_details", mode: 'copy', pattern: 'output_primer3_*'
	
	input:
	path 'in_primer3_handle' from in_primer3
	path 'primer_settings_handle' from params.primer_settings
	path 'circ_info' from circ_info
	file 'exons_bed' from exons
	val 'length' from params.temp_l
 
	output:
	path 'all_primers_circ*' into all_primers_per_circ
	path 'selected_primers_*' into results_per_circ
	path 'all_primers' into out_dir

	// this file is copied to details directory
	path 'output_primer3_*'

	"""
	/bin/primer3-2.5.0/src/primer3_core --output=output_primer3.txt --p3_settings_file=$primer_settings_handle $in_primer3_handle
	04_split_primers.py -i output_primer3.txt -l $params.design -n $length
	bedtools sort -i bed_in.txt > bed_in_sorted.txt
	bedtools closest -d -t first -a bed_in_sorted.txt -b $exons_bed > out_anno.txt
	mkdir all_primers 
	05_filter.py -A circ_info -P all_primers_circ*.txt
	"""
}

process print_output {

	publishDir params.output_dir, mode: 'copy'

	input:
	file 'results_per_circ*' from results_per_circ.collect()
	path 'all_primer_files' from out_dir.collect()

	output:
	path 'all_primers'
	path 'filtered_primers.txt'
	path 'run_summary.txt'

	"""
	mkdir all_primers
	cp all_primer_files*/* all_primers/
	echo "circ_ID	chr	start	end	strand	primer_ID	FWD_primer	FWD_rc	REV_primer	REV_rc	FWD_pos	FWD_length	FWD_pos_start	FWD_pos_end	REV_pos	REV_length	REV_pos_start	REV_pos_end	FWD_Tm	REV_Tm	FWD_GC	REV_GC	amplicon	FWD_type	REV_type	BSJ_type	same_primers	same_BSJ" > filtered_primers.txt
	cat results_per_circ* >> filtered_primers.txt
	06_make_summary.py
	"""
}
