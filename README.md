 # CIRCprimerXL_ONT
read me needs to be updated for ONT primer design!
Collaborators: Marieke Vromman, Pieter-Jan Volders

Questions concerning the GitHub structure/scripts can be addressed to any of the collaborators.

Primer design pipeline for circRNAs based on primerXL (Lefever, S., Pattyn, F., De Wilde, B. et al. High-throughput PCR assay design for targeted resequencing using primerXL. BMC Bioinformatics 18, 400 (2017). https://doi.org/10.1186/s12859-017-1809-3).

This pipeline runs entirely in the [oncornalab/primerxl_circ](https://hub.docker.com/repository/docker/oncornalab/primerxl_circ) docker image, which is available on DockerHub. It is not necessary to download this image locally, as Nextflow pulls the latest version automatically from DockerHub.

## Installation
### Running on your computer
[Nextflow](https://www.nextflow.io/) and [Docker](https://docs.docker.com/get-docker/) should be installed locally. Make sure [Docker Desktop](https://www.docker.com/products/docker-desktop) is running when you want run the pipeline.

### Running on the HPC (UGent)
Nextflow version 20.10.0 is available on all clusters (swalot, skitty, victini, joltik, kirlia, doduo). The pipeline can be run through an interactive session. The pipeline can only run from the $VSC_SCRATCH_VO_USER directory.

```
qsub -I -l nodes=1:ppn=16 -l walltime=04:00:00
cd $VSC_SCRATCH_VO_USER/CIRCprimerXL/
module load Nextflow/20.10.0
nextflow run CIRCprimerXL.nf --help
```


## General usage

```
$ nextflow run CIRCprimerXL.nf --help

Usage:
	
	The typical command for running the pipeline is as follows:
	nextflow run CIRCprimerXL.nf -profile singularity
	
	Mandatory nextflow arguments:
	-profile            set to 'local' when running locally, set to 'singularity' when running on the HPC

	Mandatory pipeline arguments:
	--input_bed         path to input file with circRNAs in bed format (0-based annotation)


	Optinal pipeline arguments:
	--splice			when set to 'yes' the input sequence will be spliced, when set to 'no' the input sequence will be unspliced
	--primer_settings   path to file with primer3plus settings (see primer3 manual)
	--chrom_file        file containing all chromosome sizes (to validate bed file)
	--primer3_diff      the minimum number of base pairs between the 3' ends of any two left primers (see also primer3 PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE)
	--primer3_nr        the number of primers designed by primer3; caution: setting this parameter to a large value will increase running time
	--min_tm	    minimum melt temperature of the primers (default: 58)
	--max_tm	    maximum melt temperature of the primers(default: 60)
	--opt_tm	    optimal melt temperature of the primers(default: 59)
	--diff_tm	    maximum difference in melt temperature between the primers(default: 2)
	--min_gc	    minimum GC contect of the primers (default: 30)
	--max_gc	    maximum GC contect of the primers(default: 80)
	--opt_gc	    optimal GC contect of the primers(default: 50)
	--amp_min	    minimum amplicon length (default: 60)
	--amp_max	    maximum amplicon length (default: 0)
	--temp_l            the number of nucleotides on each side of the circRNA BSJ that will be used for the template (example 150 => template of 300 nts in total)
	--output_dir        path to directory where the output files will be saved
	
```


You can easily create your own profiles by modifying the nextflow.config file.

Nextflow keeps track of all the processes executed in your pipeline. If you want to restart the pipeline after a bug, add '-resume'. The execution of the processes that are not changed will be skipped and the cached result used instead.

Note: The pipeline results are cached by default in the directory $PWD/work. This folder can take of lot of disk space. If your are sure you won’t resume your pipeline execution, clean this folder periodically.

Note: If a circRNA is smaller than the requested template size, the template size is reduced to the circRNA size. Of note, if this 300-nucleotide template sequence includes an exon-intron boundary, the intronic region (which may not be part of the circRNA) is included. Some circRNAs effectively also include intronic sequences, and some BSJ concatenate an exonic and an intronic sequence. 

## Output
In the output folder, you will find
<ul>
  <li>filtered_primers.txt, a file containing one selected primer pair per circRNA (see below for column details)</li>
  <li>log_file.txt, a file </li>
  <li>summary_run.txt, </li>
  <li>all_primers directory</li>
  <li>primer3_details directory</li>
</ul>

filtered_primers.txt output file column names:

| column name      | description                                                                                                            |
|:-----------------|:-----------------------------------------------------------------------------------------------------------------------|
| circ_ID          | circ id assigned to each circRNA (unique within one run)                                                               |
| chr              | circRNA chromome                                                                                                       |
| start            | circRNA start position                                                                                                 |
| end              | circRNA end position                                                                                                   |
| primer_ID        | primer ID generated by primer3                                                                                         |
| FWD_primer       | forward primer                                                                                                         |
| REV_primer       | reverse primer                                                                                                         |
| FWD_pos          | relative position of forward primer                                                                                    |
| FWD_length       | length of forward primer                                                                                               |
| REV_pos          | relative position of reverse primer                                                                                    |
| REV_length       | length of reverse primer                                                                                               |
| FWD_Tm           | melt temperature of forward primer                                                                                     |
| REV_Tm           | melt temperature of reverse primer                                                                                     |
| FWD_GC           | GC content of forward primer                                                                                           |
| REV_GC           | GC content of reverse primer                                                                                           |
| amplicon         | amplicon sequence amplified by the primer pair                                                                         |
| PASS             | result of filtering (PASS if the primer pair passed all filters, FAIL if   the primer pair failed one or more filters) |



## Other species
As default, CIRCprimerXL designs primers voor humans. To design primers for other species, the following files have to be provided and parsed through the corresponding parameters:
<ul>
  <li>a file containing the chromosome sizes (parameter chrom_file) (for example from: https://www.ncbi.nlm.nih.gov/grc/human/data)</li>
  <li>a fastahack index (parameter index_fasta) and Bowtie index (parameter index_bowtie) (see also above)</li>
  <li>a SNP database link (parameter snp_url) (for example: http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153Common.bb
)</li>
  <li>a file containing all ENST numbers of canonical transcripts or transcrits. This file can be generated by downloading the MANE file (https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/) and transforming it into a simple list using generate_ENST_list.py (can be run in the Docker image).
</ul>


## Nextflow tower

[Nextflow tower](https://tower.nf/) can be used to monitor the pipeline while it's running.
```
nextflow run CIRCprimerXL.nf -with-tower
```

When Nextflow tower is used in combination with the HPC, the nextflow version and tower access token should be indicated.
```
export NXF_VER="20.10.0"
export TOWER_ACCESS_TOKEN=your_token_here
```

