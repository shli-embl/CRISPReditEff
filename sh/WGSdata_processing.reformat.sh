#!/bin/sh
set -e
set -u
## GetOpt ##
CPU=8
usage(){
	echo "Usage: $0 [-F READ1_FILES(separated by ';')] [-R READ2_FILE(separated by ';')] [-O OUTPUT_DIR] [-S SAMPLE_NAME] [-Y YL_VERSION] [-C TARGET_CHR] [-D DONOR_START_COORD] [-d DONOR_END_COORD] [-u CPU_NUMBER (optional)]"
}
exit_abnormal(){
	usage
	exit 1
}
while getopts ":F:R:O:Y:C:D:S:d:u:" options; do
	case "${options}" in
	F)
#		SEQFILES1=${OPTARG}
		IFS=';' read -ra SEQFILES1 <<< "${OPTARG}"      ###split names (sep by ;) into array
#		N_FILE=${#SEQFILES1[@]}
		;;
	R)
#		SEQFILES2=${OPTARG}
		IFS=';' read -ra SEQFILES2 <<< "${OPTARG}" 
		;;
	O)
		OUT_DIR=${OPTARG}
		;;
	S)
		SAMPLE_NAME=${OPTARG}
		;;
	Y)
		YL=${OPTARG}
		;;
	C)
		V_CHR=${OPTARG}
		;;
	D)
		D_START=${OPTARG}
		;;
	d)
		D_END=${OPTARG}
		;;
	u)
		CPU=${OPTARG}
		;;
	:)	
		echo "Error: -${OPTARG} requires an argument."
		exit_abnormal
		;;
	*)
		exit_abnormal
		;;
	esac
done

if [[ ! -v SEQFILES1 ]] || [[ ! -v SEQFILES2 ]] || [[ ! -v OUT_DIR ]] || [[ ! -v SAMPLE_NAME ]] || [[ ! -v YL ]] || [[ ! -v V_CHR ]] || [[ ! -v D_START ]] || [[ ! -v D_END ]]; then
	echo "Error: Missing arguments."
	exit_abnormal
fi

if [[ ${#SEQFILES1[@]} -ne ${#SEQFILES2[@]} ]]; then
	echo "Error: File numbers from -F and -R are not matched."
	exit_abnormal
fi
N_FILE=${#SEQFILES1[@]}
		

#location of the reference genome and dependent in-house scripts
SH_PATH=$(dirname $0)
SRC_PATH="$SH_PATH/../src"
FASTA_PATH="$SH_PATH/../fasta"
REFERENCE_GENOME="$FASTA_PATH/yeast_reference.fa"
declare -A YL_TO_BACKBONE=(["YL151"]="pKR452" ["YL152"]="pKR452" ["YL153"]="pKR348" ["YL154"]="pKR348" ["YL155"]="pKR452" ["YL156"]="pKR452")
declare -A BACKBONE_REF_GENOME=(["pKR452"]="$FASTA_PATH/REDI_cassette_pKR452.fa" ["pKR348"]="$FASTA_PATH/REDI_cassette_pKR348.fa")
declare -A GUIDE_START_COORD=(["pKR452"]=1483 ["pKR348"]=1640)
declare -A GUIDE_END_COORD=(["pKR452"]=1502 ["pKR348"]=1659)
declare -A BC1_START_COORD=(["pKR452"]=3837 ["pKR348"]=3994)
declare -A BC1_END_COORD=(["pKR452"]=3867 ["pKR348"]=4024)
#flanking sequences of guide used for calling
declare -A GUIDE_LEFT_150BP=(["pKR452"]="AAGGTGCGCATTTTTTCACACCCTACAATGTTCTGTTCAAAAGATTTTGGTCAAACGCTGTAGAAGTGAAAGTTGGTGCGCATGTTTCGGCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATCGGAGCTGCGATTGGCAGGCGCGCC" \
["pKR348"]="TAATTTATCACTACGAAATCTTGAGATCGGGCGTTCGACTCGCCCCCGGGAGAGATGGCCGGCATGGTCCCAGCCTCCTCGCTGGCGCCGGCTGGGCAACACCTTCGGGTGGCGAATGGGACTTTGGGAGCTGCGATTGGCAGGCGCGCC")
declare -A GUIDE_RIGHT_150BP=(["pKR452"]="GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTTTTTTGTCTCGCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCA" \
["pKR348"]="GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTTTTTTGTCTCGCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCA" )
declare -A BC1_LEFT_150BP=(["pKR452"]="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCCGCCGTGGGTTTCTTCCACTCAAACATGTCATGCATCACGTGCTAGCTGACATGACGTC" \
["pKR348"]="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCCGCCGTGGGTTTCTTCCACTCAAACATGTCATGCATCACGTGCTAGCTGACATGACGTC")
declare -A BC1_RIGHT_150BP=(["pKR452"]="AGATCGGAAGAGAGTCGTGTAGGGAAAGAGCGGCCGCTGGAGGCAACAGTTACTAAACGGGTCACTTGACCGCAGTCAGCCGCATTAATGAGGGACCAGTTTAATTCGGTAATCTCCGAACAGAAGGAAGAACGAAGGAAGGAGCACAGA" \
["pKR348"]="AGATCGGAAGAGAGTCGTGTAGGGAAAGAGCGGCCGCTGGAGGCAACAGTTACTAAACGGGTCACTTGACCGCAGTCAGCCGCATTAATGAGGGACCAGTTTAATTCGGTAATCTCCGAACAGAAGGAAGAACGAAGGAAGGAGCACAGA")

# get cassette reference
BACKBONE=${YL_TO_BACKBONE[$YL]}
CASSETTE_REF_GENOME=${BACKBONE_REF_GENOME[$BACKBONE]}
# create result folders

if [ ! -d $OUT_DIR/bams ]; then
        mkdir $OUT_DIR/bams
fi

if [ ! -d $OUT_DIR/seqs ]; then
        mkdir $OUT_DIR/seqs
fi


if [ ! -d $OUT_DIR/fastqc ]; then
        mkdir $OUT_DIR/fastqc
fi


if [ ! -d $OUT_DIR/fastqc/$SAMPLE_NAME.fastqc ]; then
        mkdir $OUT_DIR/fastqc/$SAMPLE_NAME.fastqc
fi

if [ ! -d $OUT_DIR/logs ]; then
        mkdir $OUT_DIR/logs
fi

if [ ! -d $OUT_DIR/gVCF ]; then
        mkdir $OUT_DIR/gVCF
fi

if [ ! -d $OUT_DIR/VCF ]; then
        mkdir $OUT_DIR/VCF
fi

if [ ! -d $OUT_DIR/recal ]; then
        mkdir $OUT_DIR/recal
fi

if [ ! -d $OUT_DIR/map_stats ]; then
        mkdir $OUT_DIR/map_stats
fi

if [ ! -d $OUT_DIR/guide_barcode ]; then
        mkdir $OUT_DIR/guide_barcode
fi

for FILE in ${SEQFILES1[@]}
do
	if [ ! -f $FILE ]; then
        echo "Error: INPUTFILE $FILE do not exist" 
        exit_abnormal
	fi
done

for FILE in ${SEQFILES2[@]}
do
	if [ ! -f $FILE ]; then
        echo "Error: INPUTFILE $FILE do not exist" >&2
        exit_abnormal
	fi
done


#START processing reads
for ((N=0;N<N_FILE;N++))
do
	#STEP 1: TRIM ADAPTERS
	cutadapt -g AGATGTGTATAAGAGACAG -G AGATGTGTATAAGAGACAG \
		-a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
		-g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -G TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
		-a CTGTCTCTTATACACATCTGACGCTGCCGACGA -A CTGTCTCTTATACACATCTGACGCTGCCGACGA \
		-g GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -G GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
		-a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
		-o $OUT_DIR/seqs/$SAMPLE_NAME.r1.$N.trimmed.fastq -p $OUT_DIR/seqs/$SAMPLE_NAME.r2.$N.trimmed.fastq -q 20 --trim-n -m 20 --max-n=0 ${SEQFILES1[N]} ${SEQFILES2[N]}

	fastqc -o $OUT_DIR/fastqc/$SAMPLE_NAME.fastqc $OUT_DIR/seqs/$SAMPLE_NAME.r1.$N.trimmed.fastq $OUT_DIR/seqs/$SAMPLE_NAME.r2.$N.trimmed.fastq
	
	#STEP 2: MAP READS TO THE REFERENCE GENOME
	bwa mem -t $CPU $REFERENCE_GENOME -R "@RG\tID:RUN${N}\tSM:${SAMPLE_NAME}\tPL:illumina\tLB:${SAMPLE_NAME}\tPU:RUN${N}" $OUT_DIR/seqs/$SAMPLE_NAME.r1.$N.trimmed.fastq $OUT_DIR/seqs/$SAMPLE_NAME.r2.$N.trimmed.fastq | samtools view -bS > $OUT_DIR/bams/$SAMPLE_NAME.$N.bam
	
	#STEP 3: SORT BAM
	gatk SortSam -I $OUT_DIR/bams/$SAMPLE_NAME.$N.bam -O $OUT_DIR/bams/$SAMPLE_NAME.$N.sort.bam -SO coordinate --VALIDATION_STRINGENCY SILENT
	rm $OUT_DIR/bams/$SAMPLE_NAME.$N.bam
done

#STEP 4: COMBINE BAMS AND MARK PCR DUPLICATES
cat_command=""
for ((N=0;N<N_FILE;N++))
do
	cat_command="${cat_command} -I $OUT_DIR/bams/$SAMPLE_NAME.$N.sort.bam" 
done

gatk MarkDuplicates $cat_command -O $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.bam -M $OUT_DIR/bams/$SAMPLE_NAME.metrics --VALIDATION_STRINGENCY SILENT
rm $OUT_DIR/bams/$SAMPLE_NAME.*.sort.bam $OUT_DIR/bams/$SAMPLE_NAME.metrics
samtools index $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.bam

#STEP5: BASE RECALIBRATION
#RECAL1
gatk HaplotypeCaller -R $REFERENCE_GENOME -I $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.bam -ploidy 1 -O $OUT_DIR/recal/$SAMPLE_NAME.recal1.vcf
gatk BaseRecalibrator -R $REFERENCE_GENOME -I $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.bam --known-sites $OUT_DIR/recal/$SAMPLE_NAME.recal1.vcf -O $OUT_DIR/recal/${SAMPLE_NAME}.recal1.table
gatk ApplyBQSR -R $REFERENCE_GENOME -I $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.bam --bqsr-recal-file $OUT_DIR/recal/$SAMPLE_NAME.recal1.table -O $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.recal1.bam
rm $OUT_DIR/recal/$SAMPLE_NAME.recal1.vcf $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.bam  $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.bam.bai
#RECAL2
gatk HaplotypeCaller -R $REFERENCE_GENOME -I $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.recal1.bam -ploidy 1 -O $OUT_DIR/recal/$SAMPLE_NAME.recal2.vcf
gatk BaseRecalibrator -R $REFERENCE_GENOME -I $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.recal1.bam --known-sites $OUT_DIR/recal/$SAMPLE_NAME.recal2.vcf -O $OUT_DIR/recal/$SAMPLE_NAME.recal2.table
gatk ApplyBQSR -R $REFERENCE_GENOME -I $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.recal1.bam --bqsr-recal-file $OUT_DIR/recal/$SAMPLE_NAME.recal2.table -O $OUT_DIR/bams/${SAMPLE_NAME}.sort.dedup.recal2.bam
rm $OUT_DIR/recal/$SAMPLE_NAME.recal2.vcf $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.recal1.bam $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.recal1.bai $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.recal2.bai
rm $OUT_DIR/recal/$SAMPLE_NAME.recal1.vcf.idx $OUT_DIR/recal/$SAMPLE_NAME.recal2.vcf.idx
gatk AnalyzeCovariates -before $OUT_DIR/recal/$SAMPLE_NAME.recal1.table -after $OUT_DIR/recal/$SAMPLE_NAME.recal2.table -plots $OUT_DIR/recal/$SAMPLE_NAME.BQSR.pdf

#STEP6: INDEX FINAL BAM FILE; CALCULATE STATS
samtools index $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.recal2.bam
samtools idxstats $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.recal2.bam > $OUT_DIR/map_stats/$SAMPLE_NAME.idxstats
samtools flagstat $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.recal2.bam > $OUT_DIR/map_stats/$SAMPLE_NAME.flagstats

#STEP7: CLEAN THE TARGET SITE (REMOVE READS FROM LANDINGPAD)
samtools view -h $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.recal2.bam > $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.recal2.sam
perl $SRC_PATH/cleanDonorReads.pl -in $OUT_DIR/bams/$SAMPLE_NAME.sort.dedup.recal2.bam -out $OUT_DIR/bams/$SAMPLE_NAME.clean.bam -chr $V_CHR -d_start $D_START -d_end $D_END -verbose

#STEP8: CALL VARIANTS AT TARGET SITE
gatk HaplotypeCaller -R $REFERENCE_GENOME -I $OUT_DIR/bams/$SAMPLE_NAME.clean.bam -ploidy 1 -O $OUT_DIR/gVCF/$SAMPLE_NAME.vcf -ERC BP_RESOLUTION -L ${V_CHR}:${D_START}-${D_END}

gatk GenotypeGVCFs -R $REFERENCE_GENOME -V ${OUT_DIR}/gVCF/${SAMPLE_NAME}.vcf -O ${OUT_DIR}/VCF/${SAMPLE_NAME}.vcf

D_EXT_START=$((D_START - 200))
D_EXT_END=$((D_END + 200))

gatk HaplotypeCaller -R $REFERENCE_GENOME -I $OUT_DIR/bams/$SAMPLE_NAME.clean.bam -ploidy 1 -O $OUT_DIR/gVCF/$SAMPLE_NAME.ext.vcf -ERC BP_RESOLUTION -L ${V_CHR}:${D_EXT_START}-${D_EXT_END}

gatk GenotypeGVCFs -R $REFERENCE_GENOME -V ${OUT_DIR}/gVCF/${SAMPLE_NAME}.ext.vcf -O ${OUT_DIR}/VCF/${SAMPLE_NAME}.ext.vcf

#STEPS1: MAP READS TO CONSTRUCT AND CALL GUIDE AND BARCODE
for ((N=0;N<N_FILE;N++))
do
	bwa mem -t $CPU $CASSETTE_REF_GENOME $OUT_DIR/seqs/$SAMPLE_NAME.r1.$N.trimmed.fastq $OUT_DIR/seqs/$SAMPLE_NAME.r2.$N.trimmed.fastq | samtools view -bS > $OUT_DIR/bams/$SAMPLE_NAME.$N.barcode.bam
	samtools view -bF 12 $OUT_DIR/bams/$SAMPLE_NAME.$N.barcode.bam > $OUT_DIR/bams/$SAMPLE_NAME.$N.barcode.mapped.bam
rm $OUT_DIR/bams/$SAMPLE_NAME.$N.barcode.bam
	gatk SortSam -I $OUT_DIR/bams/$SAMPLE_NAME.$N.barcode.mapped.bam -O $OUT_DIR/bams/$SAMPLE_NAME.$N.barcode.sort.bam -SO coordinate --VALIDATION_STRINGENCY SILENT
rm $OUT_DIR/bams/$SAMPLE_NAME.$N.barcode.mapped.bam
done
cat_command=""
for ((N=0;N<N_FILE;N++))
do
	cat_command="${cat_command} -I $OUT_DIR/bams/$SAMPLE_NAME.$N.barcode.sort.bam" 
done

gatk MarkDuplicates $cat_command -O $OUT_DIR/bams/$SAMPLE_NAME.barcode.sort.dedup.bam -M $OUT_DIR/bams/$SAMPLE_NAME.barcode.metrics --VALIDATION_STRINGENCY SILENT
rm $OUT_DIR/bams/$SAMPLE_NAME.*.barcode.sort.bam $OUT_DIR/bams/$SAMPLE_NAME.barcode.metrics
mv $OUT_DIR/bams/$SAMPLE_NAME.barcode.sort.dedup.bam $OUT_DIR/bams/$SAMPLE_NAME.barcode.bam
samtools index $OUT_DIR/bams/$SAMPLE_NAME.barcode.bam

##calling guide and barcode
GUIDE_START=${GUIDE_START_COORD[$BACKBONE]}
GUIDE_END=${GUIDE_END_COORD[$BACKBONE]}
BC1_START=${BC1_START_COORD[$BACKBONE]}
BC1_END=${BC1_END_COORD[$BACKBONE]}
GUIDE_RIGHT=${GUIDE_RIGHT_150BP[$BACKBONE]}
GUIDE_LEFT=${GUIDE_LEFT_150BP[$BACKBONE]}
BC1_RIGHT=${BC1_RIGHT_150BP[$BACKBONE]}
BC1_LEFT=${BC1_LEFT_150BP[$BACKBONE]}
perl $SRC_PATH/extractGuideBarcode.pl -in $OUT_DIR/bams/$SAMPLE_NAME.barcode.bam -out $OUT_DIR/guide_barcode/$SAMPLE_NAME.guide.tbl -pos GRNA_BARCODEV2:$GUIDE_START-$GUIDE_END -L $GUIDE_LEFT -R $GUIDE_RIGHT 
perl $SRC_PATH/extractGuideBarcode.pl -in $OUT_DIR/bams/$SAMPLE_NAME.barcode.bam -out $OUT_DIR/guide_barcode/$SAMPLE_NAME.barcode.tbl -pos GRNA_BARCODEV2:$BC1_START-$BC1_END -L $BC1_LEFT -R $BC1_RIGHT 

for ((N=0;N<N_FILE;N++))
do
	rm $OUT_DIR/seqs/$SAMPLE_NAME.r1.$N.trimmed.fastq $OUT_DIR/seqs/$SAMPLE_NAME.r2.$N.trimmed.fastq
done
