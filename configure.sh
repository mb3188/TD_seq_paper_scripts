#!/bin/bash

###NEVER MODIFY ANYTHIG IN THE CONFIGURE FILE EXCEPT the path for FQTABLE, everythign else should be modified in the bash scripts for that specific process


#Please note that when defining paths, do not include a "/" on the right side,
#For example, define /path/to/directory but not /path/to/directory/

# ** General settings **

#note: PLOIDY=[1|2]
#      GLOBALREALIGNMENT=[0(no)|1(yes)]
#      GFF not used directly by pipeline, but will be needed for analysis scripts
#        including ko scripts and covByFeature 
#      GFF should be same GFF used to build snpEff db


#ADJUST FOR YOUR DIRECTORIES BELOW

PLOIDY=1
GLOBALREALIGNMENT=0
FQTABLE=/data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.2/fastq_to_align_MR
FQTOPDIR=/data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/fastqs
REF=/data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/BWA_genome_index/TrichDB-2.0_TvaginalisG3_Genome_cleaned.fa
GFF=/data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/BWA_genome_index/TrichDB-2.0_TvaginalisG3.gff



GENOMEVERSION=TrichDB-2.0_TvaginalisG3_Genome_cleaned


# ** callable depths **

#CALLABLEDEPTHS is an array of depth values. These depths are
#used to calculate the per sample set of intervals in the 
#aligned bam that have at least this depth of coverage.  A
#second script then calculates the intersections of these
#intervals across all samples.

#note: can be memory intensive, values should be less than 0-5 
#      for low coverage (< 15X) datasets.

CALLABLEDEPTHS=(1)


# ** set qsub -t option for parallelizing over samples **

#note:    SAMPLEWISE_T_OPT should take the form of a range 
#         from 0 to the number of genomes sequenced minus 1 
#example: if 20 genomes, then SAMPLEWISE_T_OPT=0-19 

SAMPLEWISE_T_OPT=0-15


# ** set qsub -t option for parallelizing over chromosomes **

#note:    One node is used per value (0-16=17 nodes) so please respect computational
#         limitations if you have a large fragmented genome with many contigs
#note:    CHRWISE_T_OPT take the form of a range from 0 to the number of chromosomes 
#         (contigs) minus 1

#example: if genome has 17 chromosomes, the CHRWISE_T_OPT should be 0-16

CHRWISE_T_OPT=1


## *** Settings below generally do not need to be altered *** ##

#samtools
SAMTOOLSMOD=samtools/intel/0.1.19
SAMTOOLSPATH=/share/apps/samtools/intel/0.1.19
#bwa
BWAMOD=bwa/gnu/0.7.8
BWAMODPATH=/share/apps/bwa/gnu/0.7.8

#picard-tools
PICARDMOD=picard-tools/1.111
PICARDPATH=/share/apps/picard-tools/1.111

#gatk
GATKMOD=gatk/3.1-1
GATKPATH=/share/apps/gatk/3.1-1

#pindel
PINDELMOD=pindel/intel/0.2.5a4
PINDEL2VCFMOD=/share/apps/pindel/intel/0.2.5a4

#delly
DELLYMOD=delly/intel/0.5.5

#bedtools
BEDTOOLSMOD=bedtools/intel/2.17.0

#vcftools
VCFTOOLSMOD=vcftools/intel/0.1.12a

#groupBy - when using groupby specify full path ${FILOPATH}/groupBy
FILOMOD=filo/intel/1.1.0
FILOPATH=/share/apps/filo/1.1.0/intel/bin/

#checks for some files that might be commonly mis-specified

if [ ! -f ${REF} ]; then
	echo _ERROR_ reference [ $REF ] not found
	exit 1
fi

#check that reference is indexed
if [ ! -f ${REF}.fai ]; then
        echo _ERROR_ reference does not appear to be indexed
        exit 1
fi

if [ ! -f ${FQTABLE} ]; then
        echo _ERROR_ fastq table [ ${FQTABLE} ] not found
        exit 1
fi


#validate fastq table format, validates that FASTQ table is in a proper format. 
perl -ne 'die("_ERROR_ error in fastq table format\n") unless(/^\S+\t\S+\t\S+$/)' ${FQTABLE}
