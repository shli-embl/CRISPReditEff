# CRISPReditEff
Defining CRISPR/Cas9 editing outcome based on sequencing data

Try "sh run_example.sh"!

## Commands:
### sh/WGSdata_processing.reformat.sh:
#### Description: 
Part 1) raw read processing, mapping to genome, read filtering and on-target variant calling; part 2) mapping reads to barcode locus (additional reference fasta provided), extract guide and barcode sequence. 

#### Usage: 
./WGSdata_processing.reformat.sh [-h] [-F "run1_1.fastq.gz;run2_1.fastq.gz;...;runX_1.fastq.gz] [-R "run1_2.fastq.gz;run2_2.fastq.gz;...;runX_1.fastq.gz"] [-O ouput_directory] [-S sample_name] [-Y YL_version (YL151-YL156)] [-C targeted_chromosome_id] [-D donor_start_coordinate] [-d donor_end_coordinate] [-u CPU_number(optional)]

### src/cleanDonorReads.pl:
#### Description: 
Exclude DNA donor-related reads from target sites in a given BAM file. 

#### Usage: 
perl cleanDonorReads.pl [-h] [-in input.bam] [-out output.bam] [-chr targeted_chromosome_id] [-d_start donor_start_coordinate] [-d_end donor_end_coordinate]

### src/extractGuideBarcode.pl:
#### Description: 
Identifying Guide/Barcode sequence from mapped reads (bams). The sequence on reference is represented as N strings.

#### Usage: 
perl extractGuideBarcode.pl [-h] [-in input.bam] [-out output.tbl] [-pos targeted_contig_id:start_position-end_position] [-R 3'_junction_sequence] [-L 5'_junction_sequence]
