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


echo _BEGIN_ [ extract_map_reads_nb.sh ${ACC} ]: `date`

#GENOMEINDX=$1 <- only change necessary to run script on conventional 
#desktop is to pass GENOMEINDX from command line rather than
#than set GENOMEINDX=${PBS_ARRAYID} (below)
#note: PBS_ARRAYID or GENOMEINDX should range from 0 to the number of genomes minus 1

GENOMEINDX=${PBS_ARRAYID}

GFF=/data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/BWA_genome_index/TrichDB-2.0_TvaginalisG3_gene.gff

#For trouble-shooting too many open filehandles errors:
#http://seqanswers.com/forums/showthread.php?t=7525
#http://sourceforge.net/apps/mediawiki/picard/index.php?title=Main_Page#Q:_I.27m_getting_java.io.FileNotFoundException:_.28Too_many_open_files.29._What_should_I_do.3F

#if [ ! -d ${FQTOPDIR} ]; then
	#echo _ERROR_ FQTOPDIR [ ${FQTOPDIR} ] not found
	#exit 1
#fi

#extract columns from $FQTABLE
COL1=(`cut -f 1 $FQTABLE`)
COL2=(`cut -f 2 $FQTABLE`)
COL3=(`cut -f 3 $FQTABLE`)

ACC=${COL1[$GENOMEINDX]}
FQ1=${FQTOPDIR}/${ACC}/${COL2[$GENOMEINDX]}
FQ2=${FQTOPDIR}/${ACC}/${COL3[$GENOMEINDX]}


#module purge

#modules
module load ${BWAMOD} ${SAMTOOLSMOD} ${PICARDMOD} ${GATKMOD} ${HTSEQMOD}


# this SCRIPTS extracts teh mammping statistics from each mapping file and it is run with the following command qsub -t 0-159 extract_map_reads_nb.pbs

#Htseq count for the number of reads 



samtools view -h -o /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/RNA_SEQ_ANALYSIS/project.1/alignments.1/${ACC}/bwa_msd.sam /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/RNA_SEQ_ANALYSIS/project.1/alignments.1/${ACC}/bwa_msd.bam

htseq-count --type=gene -s reverse -i "ID" /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/RNA_SEQ_ANALYSIS/project.1/alignments.1/${ACC}/bwa_msd.sam $GFF > /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/RNA_SEQ_ANALYSIS/project.1/alignments.1/${ACC}_htseq.txt

out1=`grep -c 0 /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/RNA_SEQ_ANALYSIS/project.1/alignments.1/${ACC}_htseq.txt`

htseq-count --type=gene -s reverse -a 0 -i "ID" /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/RNA_SEQ_ANALYSIS/project.1/alignments.1/${ACC}/bwa_msd.sam $GFF > /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/RNA_SEQ_ANALYSIS/project.1/alignments.1/${ACC}_htseq_no_alignment_filter.txt


out2=`grep -c 0 /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/RNA_SEQ_ANALYSIS/project.1/alignments.1/${ACC}_htseq_no_alignment_filter.txt`



#samtools flagstat ${ACC}/bwa_msd.bam > ${ACC}/bwa_msd.flagstat

#sed -n '3 p' < ${ACC}/bwa_msd.flagstat | awk '{print $1}' >> total_nb_mapped_reds.txt

#sed -n '1 p' < ${ACC}/bwa_msd.flagstat | awk '{print $1}' >> total_nb_COLLECTED_reds.txt


#awk 'FNR == 5 {print $1}' ${ACC}/bwa_msd.flagstat >> total_nb_aligned_reds.txt



out3=`awk 'FNR == 1 {print $1}' ${ACC}/bwa_msd.flagstat`
out4=`awk 'FNR == 5 {print $1}' ${ACC}/bwa_msd.flagstat` 

echo "${ACC} $out1 $out2  $out3 $out4" >> TABLE_FINAL

 

# CUFLLINKS to count the number of reads  

#cufflinks -o ${ACC}_cufflinks.txt -G $GFF ${ACC}/bwa_msd.bam 



echo _END_ [ extract_map_reads_nb.sh ${ACC} ]: `date`
