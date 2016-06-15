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
echo _BEGIN_ [ bwa aln ${ACC} ]: `date`

bwa aln -t 2 ${REF} ${FQ1} > ${ACC}/bwa_R1.sai

echo _ESTATUS_ [ R1 bwa aln ${ACC} ]: $?

bwa aln -t 2 ${REF} ${FQ2} > ${ACC}/bwa_R2.sai

echo _ESTATUS_ [ R2 bwa aln ${ACC} ]: $?

#BWA sampe
echo _BEGIN_ [ bwa sampe ${ACC} ]: `date`

bwa sampe -P -a 1300 \
-r "@RG\tID:${ACC}\tSM:${ACC}\tPL:Illumina\tLB:${ACC}\tDS:${FQ1}_${FQ2}_${DATE}" \
$REF ${ACC}/bwa_R1.sai ${ACC}/bwa_R2.sai ${FQ1} ${FQ2} \
> ${ACC}/bwa.sam

echo _ESTATUS_ [ bwa sampe ${ACC} ]: $?

#SAM->BAM conversion
samtools view -bhS ${ACC}/bwa.sam > ${ACC}/bwa.bam

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



#########################################################################################################

## I DO NOT DO PASS THIS PART BECASE I CANNOT REMOVE DUPLICATES IN DDRAD DATA BECAUSE IF I WOULD I REMOVE EVERYHTING THAT IS MORE THAN 1X BECAUSE ALL THE READS ARE STACKED ON THE SAME POSITION 
## IN THIS KIND OF SEQUNCEING BECASE OF 2 RESTRICTION ENZYMES, IF THE IS ONLY RAD SEQ DUPLicATES CAN BE REMOVED BECASE THE READS ARE NOT STACKED BECASE THERE IS SHEARING PROCESS SO 
### READS PILLE UP
########################################################################################################## 





#samtools index ${ACC}/bwa_ms.bam

#echo _ESTATUS_ [ samtools index ${ACC} bwa_ms.bam ]: $?

#ValidateSamFile
#java -Xmx28g -jar ${PICARDPATH}/ValidateSamFile.jar \
#INPUT=${ACC}/bwa_ms.bam \
#OUTPUT=${ACC}/bwa_ms.validate \
#MODE=SUMMARY \
#TMP_DIR=${TMPDIRPATH} \
#MAX_OPEN_TEMP_FILES=2000 \
#MAX_OUTPUT=1000000000 \
#MAX_RECORDS_IN_RAM=10000000

#echo _ESTATUS_ [ ValidateSamFile bwa_ms.bam ${ACC} ]: $?





#echo _BEGIN_ [ MarkDuplicates ${ACC} ]: `date`

#MarkDuplcates
#note: java.io.FileNotFoundException (Too many open files) error
#      was a problem until I implemented some changes to 
#      ValidateSam and MarkDups.  Note that on 2.3.13 I reduced 
#      MAX_FILE_HANDLE_FOR_READ_ENDS_MAP in MarkDups 
#      to 1000. I also reduced the MAX_OPEN_TEMP_FILES in 
#      ValidateSam to 1000, and increased MAX_RECORDS_IN_RAM to 
#      20,000,000 in both. This is in addition to always 
#      specifying TMP_DIR as a directory in scratch.

#note: marking duplicates returns a class of unmapped reads
#      whose mates are duplicates and mapped. When removing
#      duplicates this leads to the unmapped reads being
#      file orphans.  We remove this class of unmapped reads
#      while retaining unmapped reads not paired with 
#      duplicates (which are used in downstream split-read
#      analysis).  See procedure below.

#note: mark the duplicate reads (reads not removed at this step
#      as of pipeline version 2.3)

#java -Xmx28g -jar ${PICARDPATH}/MarkDuplicates.jar \
#INPUT=${ACC}/bwa_ms.bam \
#OUTPUT=${ACC}/bwa_msd_marked.bam \
#METRICS_FILE=${ACC}/bwa_msd_marked.metrics \
#REMOVE_DUPLICATES=false \
#ASSUME_SORTED=true \
#VALIDATION_STRINGENCY=LENIENT \
#READ_NAME_REGEX="[-a-zA-Z0-9]+\:[a-zA-Z0-9]+\:[a-zA-Z0-9]+\:[0-9]+\:([0-9]+)\:([0-9]+)\:([0-9]+)" \
#TMP_DIR=${TMPDIRPATH} \
#MAX_RECORDS_IN_RAM=1000000 \
#MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

#echo _ESTATUS_ [ MarkDuplicates ${ACC} ]: $?

#Failure to produce bwa_msd_marked.bam is common hpc failure.
#Die here if .bam not produced

#if [ ! -f ${ACC}/bwa_msd_marked.bam ]; then
	#echo _ERROR_ bwa_msd_marked.bam not produced for ${ACC}
	#exit 1
#fi

#cleanup
#rm ${ACC}/bwa_ms.bam ${ACC}/bwa_ms.bam.bai

#samtools index ${ACC}/bwa_msd_marked.bam

#echo _ESTATUS_ [ samtools index bwa_msd_marked.bam ${ACC} ]: $?

#get counts
#note: dups with mate mapped cannot be calculated directly from
#      SAM flag

#TOTREADS=`samtools view -c ${ACC}/bwa_msd_marked.bam`
#TOTDUPS=`samtools view -c -f 1024 ${ACC}/bwa_msd_marked.bam`
#DUPSMATEUNMAPPED=`samtools view -c -f 1032 ${ACC}/bwa_msd_marked.bam`
#UNMAPPEDDUPS=`samtools view -c -f 1028 ${ACC}/bwa_msd_marked.bam`

#echo _COUNTS_ Total Reads [ ${ACC} ]: ${TOTREADS}
#echo _COUNTS_ Total Duplicate Reads [ ${ACC} ]: ${TOTDUPS}
#echo _COUNTS_ Duplicate reads with mates unmapped [ ${ACC} ]: ${DUPSMATEUNMAPPED}
#echo _COUNTS_ Unmapped reads that are duplicates [ ${ACC} ]: ${UNMAPPEDDUPS}

#make list of duplicate reads whose mates are unmapped
#samtools view -f 1032 ${ACC}/bwa_msd_marked.bam | cut -f1 > ${ACC}/bwa_msd_marked.f1032.readsToExclude

#echo _ESTATUS_ [ samtools view 1032 readsToExclude ${ACC} ]: $?




#echo _BEGIN_ [ remove duplicates ${ACC} ]: `date`

#remove duplicates
#samtools view -F 1024 -bh ${ACC}/bwa_msd_marked.bam > ${ACC}/bwa_msd_marked_duprm.bam

#echo _ESTATUS_ [ samtools view -F 1024 in ${ACC} ]: $?

#cleanup
#rm ${ACC}/bwa_msd_marked.bam ${ACC}/bwa_msd_marked.bam.bai

#samtools index ${ACC}/bwa_msd_marked_duprm.bam

#echo _ESTATUS_ [ samtools index bwa_msd_marked_duprm.bam {$ACC} ]: $?

#TOTREADSAFTERDUPRM=`samtools view -c ${ACC}/bwa_msd_marked_duprm.bam`

#echo _COUNTS_ Total reads after samtools view -F 1024 [ ${ACC} ]: ${TOTREADSAFTERDUPRM}




#echo _BEGIN_ [ FilterSamReads ${ACC} ]: `date`

#remove unmapped reads whose mates are duplicated using .readsToExclude file

#java -Xmx28g -jar ${PICARDPATH}/FilterSamReads.jar \
#INPUT=${ACC}/bwa_msd_marked_duprm.bam \
#EXCLUDE_READS=${ACC}/bwa_msd_marked.f1032.readsToExclude \
#WRITE_READS_FILES=false \
#SORT_ORDER=coordinate \
#OUTPUT=${ACC}/bwa_msd.bam \
#TMP_DIR=${TMPDIRPATH} \
#VALIDATION_STRINGENCY=LENIENT \
#MAX_RECORDS_IN_RAM=1000000

#echo _ESTATUS_ [ FilterSamReads ${ACC} ]: $?

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
#rm ${ACC}/bwa_msd.bam ${ACC}/bwa_msd.bam.bai

#remove GATK-formatted index which is automatically produced by IndelRealigner
#rm ${ACC}/bwa_msdr.bai

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
samtools flagstat ${ACC}/bwa_msdr.bam > ${ACC}/bwa_msdr.flagstat

echo _ESTATUS_ [ samtools flagstat ${ACC} ]: $?

#create .md5 file with md5 checksum
md5sum ${ACC}/bwa_msdr.bam > ${ACC}/bwa_msdr.bam.md5

#create bam with only unmapped reads (for de novo assembly)
samtools view -bh -f 4 ${ACC}/bwa_msdr.bam > ${ACC}/bwa_msdr_unmapped.bam

echo _ESTATUS_ [ bwa_msdr_unmapped.bam samtools view ${ACC} ]: $?

#cleanup
#note: .sai files are large so remove bc never used downstream

#rm ${ACC}/bwa_R1.sai
#rm ${ACC}/bwa_R2.sai
rm -r ${TMPDIRPATH}

echo _END_ [ aln_pipeline.sh ${ACC} ]: `date`
