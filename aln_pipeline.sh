#!/bin/bash

#define pipeline shell variables
if [ -f ../configure.sh ]; then
	source ../configure.sh
else
	echo _ERROR_ configure.sh not found
	exit 1
fi


#Description: This script aligns paired-end fastqs to a 
#             reference genome and conducts basic
#             post-processing of alignments.

#note:        Reads are not filtered (removed) based on MAPQ
#             Optical and PCR duplicates are removed


echo _BEGIN_ [ aln_pipeline.sh ${ACC} ]: `date`

#GENOMEINDX=$1 <- only change necessary to run script on conventional 
#desktop is to pass GENOMEINDX from command line rather than
#than set GENOMEINDX=${PBS_ARRAYID} (below)
#note: PBS_ARRAYID or GENOMEINDX should range from 0 to the number of genomes minus 1

GENOMEINDX=${PBS_ARRAYID}

#For trouble-shooting too many open filehandles errors:
#http://seqanswers.com/forums/showthread.php?t=7525
#http://sourceforge.net/apps/mediawiki/picard/index.php?title=Main_Page#Q:_I.27m_getting_java.io.FileNotFoundException:_.28Too_many_open_files.29._What_should_I_do.3F

if [ ! -d ${FQTOPDIR} ]; then
	echo _ERROR_ FQTOPDIR [ ${FQTOPDIR} ] not found
	exit 1
fi

#extract columns from $FQTABLE
COL1=(`cut -f 1 $FQTABLE`)
COL2=(`cut -f 2 $FQTABLE`)
COL3=(`cut -f 3 $FQTABLE`)

ACC=${COL1[$GENOMEINDX]}
FQ1=${FQTOPDIR}/${ACC}/${COL2[$GENOMEINDX]}
FQ2=${FQTOPDIR}/${ACC}/${COL3[$GENOMEINDX]}

#make a directory for each accession listed in column 1 of $FQTABLE
if [ -d ${ACC} ]; then
	echo _ERROR_ directory ${ACC} already exists
	exit 1
else
	mkdir ${ACC}
fi

if [ ! -f ${FQ1} ]; then
	echo _ERROR_ FQ1 [ ${FQ1} ] not found
	exit 1
fi

if [ ! -f ${FQ2} ]; then
	echo _ERROR_ FQ2 [ ${FQ2} ] not found
	exit 1
fi


MKTEMPVERSION=`mktemp -V | head -1`
if [ "${MKTEMPVERSION}" != "mktemp (GNU coreutils) 8.4" ]; then
        echo _ERROR_ check compatibility of mktemp version
        exit 1
else
        TMPDIRPATH=`mktemp -d -p ${PWD}` || exit 1
fi


#modules
module load ${BWAMOD} ${SAMTOOLSMOD} ${PICARDMOD} ${GATKMOD}

echo `module list`
echo `set`
echo ulimit: `ulimit -n`
echo file descriptors: `cat /proc/sys/fs/file-max`
echo Sample: ${ACC}
echo Temporary directory: ${TMPDIRPATH}

#BWA aln with -t number of threads - be sure at least this many are requested in pbs submission script


#First do the individual read mapping with bwa bwasw

echo _BEGIN_ [ bwa bwasw ${ACC} ]: `date`


bwa bwasw -t 2 ${REF} ${FQ1} > ${ACC}/bwa_R1.sam

echo _ESTATUS_ [ R1 bwa bwasw ${ACC} ]: $?

bwa bwasw -t 2 ${REF} ${FQ2} > ${ACC}/bwa_R2.sam

echo _ESTATUS_ [ R2 bwa bwasw ${ACC} ]: $?

#BWA mem paired end
echo _BEGIN_ [ bwa mem ${ACC} ]: `date`

bwa mem -t 2 \
-aM -I 500,500,900,200 -R "@RG\tID:${ACC}\tSM:${ACC}\tPL:Illumina\tLB:${ACC}\tDS:${FQ1}_${FQ2}_${DATE}" \
$REF  ${FQ1} ${FQ2} \
> ${ACC}/bwa.sam

echo _ESTATUS_ [ bwa mem ${ACC} ]: $?

#SAM->BAM conversion

samtools view -bhS ${ACC}/bwa.sam > ${ACC}/bwa.bam

samtools view -bhS ${ACC}/bwa_R1.sam > ${ACC}/bwa_R1.bam

samtools view -bhS ${ACC}/bwa_R2.sam > ${ACC}/bwa_R2.bam

echo _ESTATUS_ [ sam to bam samtools view ${ACC} ]: $?


#cleanup SAM
#rm ${ACC}/bwa.sam

#note: we do not index here bc .bam is not coordinate-sorted


echo _BEGIN_ [ FixMateInformation ${ACC} ]: `date`

#FixMateInformation
java -Xmx28g -jar ${PICARDPATH}/FixMateInformation.jar \
INPUT=${ACC}/bwa.bam \
OUTPUT=${ACC}/bwa_m.bam \
TMP_DIR=${TMPDIRPATH} \
VALIDATION_STRINGENCY=LENIENT \
MAX_RECORDS_IN_RAM=10000000

echo _ESTATUS_ [ FixMateInformation ${ACC} ]: $?

#cleanup - no .bai file for bwa.bam to remove
#rm ${ACC}/bwa.bam

#note: we do not index here bc .bam is not coordinate-sorted

echo _BEGIN_ [ ValidateSamFile bwa_m.bam ${ACC} ]: `date`

#ValidateSamFile
java -Xmx28g -jar ${PICARDPATH}/ValidateSamFile.jar \
INPUT=${ACC}/bwa_m.bam \
OUTPUT=${ACC}/bwa_m.validate \
MODE=SUMMARY \
TMP_DIR=${TMPDIRPATH} \
MAX_OPEN_TEMP_FILES=2000 \
MAX_OUTPUT=1000000000 \
MAX_RECORDS_IN_RAM=10000000

echo _ESTATUS_ [ ValidateSamFile bwa_m.bam ${ACC} ]: $?


echo _BEGIN_ [ SortSam ${ACC} ]: `date`

#SortSam
#note: IGNORE option won't work with SamSort
#(eg IGNORE=INVALID_MAPPING_QUALITY \ IGNORE=INVALID_CIGAR) so use VALIDATION_STRINGENCY=LENIENT

java -Xmx28g -jar ${PICARDPATH}/SortSam.jar \
INPUT=${ACC}/bwa_m.bam \
OUTPUT=${ACC}/bwa_msd.bam \
SORT_ORDER=coordinate \
TMP_DIR=${TMPDIRPATH} \
MAX_RECORDS_IN_RAM=10000000 \
VALIDATION_STRINGENCY=LENIENT

echo _ESTATUS_ [ SortSam ${ACC} ]: $?

#cleanup - no .bai file for bwa_m.bam to remove
#rm ${ACC}/bwa_m.bam



#cleanup
#rm ${ACC}/bwa_msd_marked_duprm.bam ${ACC}/bwa_msd_marked_duprm.bam.bai

samtools index ${ACC}/bwa_msd.bam

echo _ESTATUS_ [ samtools index bwa_msd.bam ${ACC} ]: $?

TOTREADSAFTERFILTERSAM=`samtools view -c ${ACC}/bwa_msd.bam`

echo _COUNTS_ Total reads after FilterSamReads [ ${ACC} ]: ${TOTREADSAFTERFILTERSAM} 

echo _BEGIN_ [ ValidateSamFile bwa_msd.bam ${ACC} ]: `date`

#ValidateSamFile
java -Xmx28g -jar ${PICARDPATH}/ValidateSamFile.jar \
INPUT=${ACC}/bwa_msd.bam \
OUTPUT=${ACC}/bwa_msd.validate \
MODE=SUMMARY \
TMP_DIR=${TMPDIRPATH} \
MAX_OPEN_TEMP_FILES=2000 \
MAX_OUTPUT=1000000000 \
MAX_RECORDS_IN_RAM=10000000

echo _ESTATUS_ [ ValidateSamFile bwa_msd.bam ${ACC} ]: $?

#making assumption that markduplicates will output a sorted BAM
#given that its input was sorted.
#**can check this by dumping header of bwa_msd.bam and checking
#SO tag in header (since Picard suite honors the standard of 
#setting this tag)



echo _BEGIN_ [ RealignerTargetCreator ${ACC} ]: `date`

#RealignerTargetCreator
java -Xmx26g -jar ${GATKPATH}/GenomeAnalysisTK.jar \
-I ${ACC}/bwa_msd.bam \
-R ${REF} \
-T RealignerTargetCreator \
-o ${ACC}/bwa_msd.bam.intervals

echo _ESTATUS_ [ RealignerTargetCreator ${ACC} ]: $?

echo _BEGIN_ [ IndelRealigner ${ACC} ]: `date`

#IndelRealigner - automatically creates *.bai index for output BAM
java -Xmx26g -jar ${GATKPATH}/GenomeAnalysisTK.jar \
-I ${ACC}/bwa_msd.bam \
-R ${REF} \
-T IndelRealigner \
-targetIntervals ${ACC}/bwa_msd.bam.intervals \
-o ${ACC}/bwa_msdr.bam

echo _ESTATUS_ [ IndelRealigner ${ACC} ]: $?

#note: do not need to run fix-mate information again after re-alignment.   
#      see:  http://gatkforums.broadinstitute.org/discussion/1562/need-to-run-a-step-with-fixmateinformation-after-realignment-step

#note: re-sorting BAM is not necessary after IndelRealigner
#      see: http://gatkforums.broadinstitute.org/discussion/comment/5151#Comment_5151

#cleanup
rm ${ACC}/bwa_msd.bam ${ACC}/bwa_msd.bam.bai

#remove GATK-formatted index which is automatically produced by IndelRealigner
rm ${ACC}/bwa_msdr.bai

#Create samtools bam index
samtools index ${ACC}/bwa_msdr.bam

echo _ESTATUS_ [ samtools index bwa_msdr.bam ${ACC} ]: $?

echo _BEGIN_ [ ValidateSamFile bwa_msdr.bam ${ACC} ]: `date`

#ValidateSamFile
java -Xmx18g -jar ${PICARDPATH}/ValidateSamFile.jar \
INPUT=${ACC}/bwa_msdr.bam \
OUTPUT=${ACC}/bwa_msdr.validate \
MODE=SUMMARY \
TMP_DIR=${TMPDIRPATH} \
MAX_OPEN_TEMP_FILES=4000 \
MAX_OUTPUT=1000000000 \
MAX_RECORDS_IN_RAM=10000000

echo _ESTATUS_ [ ValidateSamFile bwa_msdr.bam ${ACC} ]: $?




#create flagstat
samtools flagstat ${ACC}/bwa.bam > ${ACC}/bwa_PAIRED.flagstat

samtools flagstat ${ACC}/bwa_R1.bam > ${ACC}/bwa_R1.flagstat

samtools flagstat ${ACC}/bwa_R2.bam > ${ACC}/bwa_R2.flagstat

echo _ESTATUS_ [ samtools flagstat ${ACC} ]: $?

#create .md5 file with md5 checksum
md5sum ${ACC}/bwa_msdr.bam > ${ACC}/bwa_msdr.bam.md5

#create bam with only unmapped reads (for de novo assembly)
samtools view -bh -f 4 ${ACC}/bwa_msdr.bam > ${ACC}/bwa_msdr_unmapped.bam

echo _ESTATUS_ [ bwa_msdr_unmapped.bam samtools view ${ACC} ]: $?

#cleanup
#note: .sai files are large so remove bc never used downstream

rm ${ACC}/bwa_R1.sai
rm ${ACC}/bwa_R2.sai
rm -r ${TMPDIRPATH}

echo _END_ [ aln_pipeline.sh ${ACC} ]: `date`
