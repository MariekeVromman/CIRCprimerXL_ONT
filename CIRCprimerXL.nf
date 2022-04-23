#!/usr/bin/env nextflow

// set default parameters

params.primer_settings = "$baseDir/assets/primer3plus_settings.txt"
params.chrom_file = "$baseDir/assets/GRCh38/chrom_sizes_GRCh38.txt"
params.input_bed = "example/path"
params.known_exons = "$baseDir/assets/GRCh38/known_exons_GRCh38.bed"

params.primer3_diff = 1
params.primer3_nr = 20
params.min_tm = 58
params.max_tm = 60
params.opt_tm = 59
params.diff_tm = 2
params.min_gc = 30
params.max_gc = 80
params.opt_gc = 50
params.amp_min = 50
params.amp_max = 0 // this param is set to 0, so that it can be adjusted depending on temp_l if the user does not supply amp_max
params.temp_l = 150


// required parameters
input_bed = file(params.input_bed)
chrom_file = file(params.chrom_file)

// print run info

log.info """\
==============================================
CIRCprimerXL pipeline ONT
==============================================
OncoRNALab - Marieke Vromman
https://github.com/OncoRNALab/CIRCprimerXL
https://hub.docker.com/repository/docker/oncornalab/CIRCprimerXL
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
	path 'all_circ.txt' into all_cic

	"""
	validate_bed.py -i $input_bed_handle -c $chrom_file_handle
	split_circRNAs.py -i $input_bed_handle
	"""

}


process get_seq {
	input:
	file ind_circ_file_handle from ind_circ_file.flatten()
	val 'diff' from params.primer3_diff
	val 'nr' from params.primer3_nr
	val 'length' from params.temp_l

	output:
	path 'input_primer3*' into in_primer3
	path 'input_filter*' into in_filter

	"""
	get_circ_seq.py -n $length -i $ind_circ_file_handle -d $diff -p $nr -m 'marieke.vromman@ugent.be'
	"""
}



process get_primers {

	publishDir "$params.output_dir/primer3_details", mode: 'copy', pattern: 'output_primer3_*'
	
	input:
	path 'in_primer3_handle' from in_primer3
	val 'temp_l_handle' from params.temp_l
	path 'primer_settings_handle' from params.primer_settings
	path 'in_filter_handle' from in_filter
 
	output:
	path 'all_primers_circ*' into all_primers_per_circ
	path 'output_primer3*'
	path 'selected_primers_*' into results_per_circ
	path('all_primers') into out_dir

	"""
	/bin/primer3-2.5.0/src/primer3_core --output=output_primer3.txt --p3_settings_file=$primer_settings_handle $in_primer3_handle
	split_primers.py -i output_primer3.txt
	mkdir all_primers
	filter.py -A in_filter_handle -P all_primers_circ*.txt -l $temp_l_handle
	gather_output.py -i all_primers/filtered_primers_*
	"""
}

process print_output {

	publishDir params.output_dir, mode: 'copy'

	input:
	file 'results_per_circ*' from results_per_circ.collect()
	path 'all_primer_files' from out_dir.collect()
	path 'all_circ_file' from all_cic


	output:
	path 'all_primers'
	path 'filtered_primers.txt'

	"""
	mkdir all_primers
	cp all_primer_files*/* all_primers/
	echo "circ_ID	chr	start	end	primer_ID	FWD_primer	FWD_rc	REV_primer	REV_rc	FWD_pos	FWD_length	REV_pos	REV_length	FWD_Tm	REV_Tm	FWD_GC	REV_GC	amplicon	PASS	start_annotation	end_annotation	splicing" > filtered_primers.txt
	cat results_per_circ* >> filtered_primers.txt
	"""
}
