#!/bin/bash
#echo _BEGIN_ [ test.sh ]: `date`



GENOMEINDX=${PBS_ARRAYID}


FQTABLE=/data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.2/fastq_to_align_MR

COL1=(`cut -f 1 $FQTABLE`) #
COL2=(`cut -f 2 $FQTABLE`) # 
COL3=(`cut -f 3 $FQTABLE`) #


ACC=${COL1[$GENOMEINDX]}

FQ1=${COL2[$GENOMEINDX]}
FQ2=${COL3[$GENOMEINDX]}



#+++++++++++++++++++++++++++++++++++++++++++++  FILTERING HIgh QUALITY READS ADN RELIABLE MAPPING FOR EACH STRAIN 

#IN order to get relible alignment I filter first read two R2.bam that have TE sequnce for mapping qualit of 30 and more

module load samtools
module load bedtools

# READs quility filtering 
#samtools view -q 30 -b /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.1/alignments.1/${ACC}/bwa_R2.bam > /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.1/alignments.1/alignments.1/${ACC}/bwa_R2.Q30.bam

# REMOVE irrelevant files

rm /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.2/alignments.1/${ACC}/bwa_R2.Q30.10x_NEW_INSERTIONS
rm /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.2/alignments.1/${ACC}/bwa_R2.Q30.5x_NEW_INSERTIONS

rm /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.2/alignments.1/${ACC}/bwa_R2.Q30.5x_1000bp_adjusted.bed

#Sort bam file
#samtools sort /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.1/alignments.1/${ACC}/bwa_R2.Q30.bam /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.1/alignments.1/${ACC}/bwa_R2.Q30_sorted



# Next I filter it also fo depth, filter for read depth of 5x and 10x, so will decide later what will be used
#bedtools genomecov -ibam /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.1/alignments.1/${ACC}/bwa_R2.Q30_sorted.bam -bg| awk '$4>10' | awk '{print $1, $2,$3}' | tr ' ' \\t > /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.1/alignments.1/${ACC}/bwa_R2.Q30.10x.bed

#bedtools genomecov -ibam /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.1/${ACC}/bwa_R2.Q30_sorted.bam -bg| awk '$4>5' | awk '{print $1, $2,$3}' | tr ' ' \\t > /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.1/alignments.1/${ACC}/bwa_R2.Q30.5x.bed


# Sort the BED file that will be used for TE INSERTION DETECTION

sort -k1,1 -k2,2n /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.2/alignments.1/${ACC}/bwa_R2.Q30.10x.bed > /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.2/alignments.1/${ACC}/bwa_R2.Q30.10x.sorted.bed

sort -k1,1 -k2,2n /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.2/alignments.1/${ACC}/bwa_R2.Q30.10x.bed > /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.2/alignments.1/${ACC}.bwa_R2.Q30.10x.sorted.bed


bedtools merge -d 300 -i /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.2/alignments.1/${ACC}.bwa_R2.Q30.10x.sorted.bed > /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.2/alignments.1/${ACC}.bwa_R2.Q30.10x.sorted.mereged.bed

#++++++++++++++++++++++++++++++++++++++++++++ MAKING GFF BED FILE 

#Now I have bed file which tell me the intervals where the reads map with qulity of 30 and 10x or 5x
# I use this bed file now to intersect with GFF.bed sorted file that I have created in order to print out all the geens that first upstream or downstream to where the insertion occured

#make GFF. bed file
# preparation of GFF files to work with bedtools closest: GFF has to be transformed into bed file
#GFF file needs to have genes only and format should be CONTIG_name  Position1   Postition2  GENE_INFO
#After that GFF.bed has to be sorted and than I can perform evaluation with bedtools closest

#make gff.bed file, replace = and ; in the gff file so you can parse the file


#sed -i -e 's/=/ /g' TrichDB-2.0_TvaginalisG3_gene.gff

#sed -i -e 's/;/ /g' TrichDB-2.0_TvaginalisG3_gene.gff


#Extracting only columns that are important, also extract the strand

#awk '{print $1, $4,$5,$7, $12, $14}' TrichDB-2.0_TvaginalisG3_gene.gff | tr ' ' \\t > TrichDB-2.0_TvaginalisG3_gene.GFF.bed

# Sort the GFF.bef file

#sort -k1,1 -k2,2n TrichDB-2.0_TvaginalisG3_gene.GFF.bed >TrichDB-2.0_TvaginalisG3_gene.GFF.sorted.bed



#++++++++++++++++++++++++++++++++++++++++++++++    FINDING WHERE ALL INDERTIONS ARE IN RELATIONSHIP TO A GENES IN THE GENOME
# using -D a which means Report distance with respect to A. When A is on the - strand, “upstream” means B has a higher (start,stop)
# as described here: http://bedtools.readthedocs.org/en/latest/content/tools/closest.html

#bedtools closest -a /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.1/TrichDB-2.0_TvaginalisG3_gene.GFF.sorted.bed -b /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.2/alignments.1/${ACC}/bwa_R2.Q30.10x.sorted.bed -D a > /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.2/alignments.1/${ACC}.10x.insertions


# If you want ot run across all isolates here is how you do it but not with this script, you should run it on the commnad line


# Extract insertions in Mariner

#bedtools closest -a /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.1/TrichDB-2.0_TvaginalisG3_gene.GFF.sorted.bed -b B7268_20a_MR.bwa_R2.Q30.10x.sorted.bed B7268_metR_MR.bwa_R2.Q30.10x.sorted.bed BUSH04_MR.bwa_R2.Q30.10x.sorted.bed BUSH20_MR.bwa_R2.Q30.10x.sorted.bed C1_MR.bwa_R2.Q30.10x.sorted.bed CHAR37_MR.bwa_R2.Q30.10x.sorted.bed G3_MR.bwa_R2.Q30.10x.sorted.bed GOR69_MR.bwa_R2.Q30.10x.sorted.bed JAM20_MR.bwa_R2.Q30.10x.sorted.bed MOR31_MR.bwa_R2.Q30.10x.sorted.bed PMGH_25_MR.bwa_R2.Q30.10x.sorted.bed SAA71_16_MR.bwa_R2.Q30.10x.sorted.bed SD2_MR.bwa_R2.Q30.10x.sorted.bed T1_MR.bwa_R2.Q30.10x.sorted.bed TAIHS176_MR.bwa_R2.Q30.10x.sorted.bed -names B7268_20 B7268_metR BUSH04 BUSH20 C1 CHAR37 G3 GOR69 JAM20 MOR31 PMGH_25 SAA71_16 SD2 T1 TAIHS176 -D a > TABLE_each_ISOLATE_reported

#bedtools closest -a /data/cgsb/carlton/MB_ddRAD/MB_BACKUP/TD_SEQ_PAPAER_ANALYSIS/TD_SEQ_ANALYSIS/project.1/TrichDB-2.0_TvaginalisG3_gene.GFF.sorted.bed -b B7268_20a_MR.bwa_R2.Q30.10x.sorted.bed B7268_metR_MR.bwa_R2.Q30.10x.sorted.bed BUSH04_MR.bwa_R2.Q30.10x.sorted.bed BUSH20_MR.bwa_R2.Q30.10x.sorted.bed C1_MR.bwa_R2.Q30.10x.sorted.bed CHAR37_MR.bwa_R2.Q30.10x.sorted.bed G3_MR.bwa_R2.Q30.10x.sorted.bed GOR69_MR.bwa_R2.Q30.10x.sorted.bed JAM20_MR.bwa_R2.Q30.10x.sorted.bed MOR31_MR.bwa_R2.Q30.10x.sorted.bed PMGH_25_MR.bwa_R2.Q30.10x.sorted.bed SAA71_16_MR.bwa_R2.Q30.10x.sorted.bed SD2_MR.bwa_R2.Q30.10x.sorted.bed T1_MR.bwa_R2.Q30.10x.sorted.bed TAIHS176_MR.bwa_R2.Q30.10x.sorted.bed -names B7268_20 B7268_metR BUSH04 BUSH20 C1 CHAR37 G3 GOR69 JAM20 MOR31 PMGH_25 SAA71_16 SD2 T1 TAIHS176 -mdb all -D a > TABLE_all_ISOLATE_reported


exit


