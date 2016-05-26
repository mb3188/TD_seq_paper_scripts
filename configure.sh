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


PLOIDY=1
GLOBALREALIGNMENT=0
FQTABLE=/data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.2/fastq_to_align_MR
FQTOPDIR=/data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/fastqs
REF=/data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/BWA_genome_index/TrichDB-2.0_TvaginalisG3_Genome_cleaned.fa
GFF=/data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/BWA_genome_index/TrichDB-2.0_TvaginalisG3.gff



# ** Genome Version **

#GENOMEVERSION is used by snpEff.  It should be the basename of your
#fasta-formatted reference (e.g., basename.fa), and should be 
#the same as (1) the "genome name" provided in the "genome:" section
#your v2.0.5 snpEff.config file, and (2) the basename of your 
#basename.genome entry in the "Genome details" section of your 
#v2.0.5 snpEff.config     (see also the pipeline README)

#Example: for GENOMEVERSION=Creinhardtii_236_cm, the 
#reference would be Creinhardtii_236_cm.fa, the "genome:" entry would
#in snpEff.config be Creinhardtii_236_cm, and the "Genome details"
# entry would be Creinhardtii_236_cm.genome


#I should crate new genome version for SNPEFF, as now Tvag 2.0 genome version is released. 

GENOMEVERSION=TrichDB-2.0_TvaginalisG3_Genome_cleaned


# ** snpEff **

#note: set path to snpEff directory (which contains the snpEff.jar) and 
#snpEff.config these will normally be in your home directory

SNPEFFPATH=/home/mb3188/snpEff
SNPEFFCONFIG=/home/mb3188/snpEff/snpEff.config


# ** snp database for base quality recalibration **

#note: only used in base recalibration step (bqsr_master.sh). May
#      leave undefined if snpdb is created by an initial round of
#      of snp-calling using the pipeline (i.e., when there is
#      not a snpdb available).
#note: SNPDB will usually be a full path to a vcf file. Please 
#      check GATK BaseRecalibrator documentation if your snpdb
#      is in a different format.

SNPDB=


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

#snpeff
#note: no backward comptability with databases when install new version, so need to remake when upgrade to 3.0
#note: if upgrade to v3.0, need to change -o vcf to -o gatk in in snpEff_pipeline.sh (see snpEff FAQ)
#note: if use new genome, $GENOME must match the name in the snpEff data dir
#note: notice snpEff.config has data_dir var that provides path to snpeff data dir.
#         may also need to be adjusted when using new genome
#note: snpEff_pipeline.sh uses -treatAllAsProteinCoding Auto. This will cause
#         all transcripts to be treated as being protein coding (true for chlamy and datepalm?) 
#         May need to be adjusted for other genomes.

#Set the genome version for snpEff.  This is the (geneome).fa : something in the snpEff.config file
#From the snpEff.config file:
## Databases & Genomes
#
# One entry per genome version.
#
# For genome version 'ZZZ' the entries look like
#	ZZZ.genome              : Real name for ZZZ (e.g. 'Human')

#  Example: In the "genome" section of the snpEff.config
#  genome: \
#  , GRCh37.65 \  #<- this is an example entry

#  Same example, in the section beginning:
#  #---
#  # Genome details
#  #---
#  #ENSEMBL v65 : Homo_sapiens.GRCh37.65
#  GRCh37.65.genome : Homo_sapiens

#The snpEff GENOMEVERSION variable would be GRCh37.65 for the above example.

#note: don't upgrade snpeff unless first verify that GATK will take new version of snpeff output (which it wouldn't as of 2.10.13)





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

if [ ! -f ${SNPEFFCONFIG} ]; then
	echo _ERROR_ snpEff.config file [ ${SNPEFFCONFIG} ] not found
	exit 1
fi

if [ ! -f ${SNPEFFPATH}/snpEff.jar ]; then
	echo _ERROR_ snpEff.jar [ ${SNPEFFPATH}/snpEff.jar ] not found
	exit 1
fi

#validate fastq table format
#perl -ne 'die("_ERROR_ error in fastq table format\n") unless(/^\S+\t\S+\t\S+$/)' ${FQTABLE}
