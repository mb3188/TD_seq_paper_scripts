** Introduction **

These are the quick and dirty instructions.  

The  pipeline provides a semi-automated method for Transposable element (TE) insertion 
variant-calling from paired-end Miseq 2*350bp data with a reference genome available. 
First step of the pipeline aligns reads to reference genome. Second phase filteres reads based
on quality control.  

The pipeline is run in two phases (1) the initial alignment for long reads
(2) a second phase where quality control, coverage filtering and insertion detection.   
The two phases are executed separately by running two scripts -- master.sh  
which produces output in alignments.1 directory and  INSERTION_DETECTION.pbs which should be copied to 
 alignments.1 directory together with its bash script INSERTION_DETECTION.pb and run after alignment has finished. 

master.sh conducts alignments 

The recommended protocol is to run aln_pipeline.pbs to completion and then execute
INSERTION_DETECTION.pbs (see below). 


 
Both master.sh and bqsr_master.sh are configured by modifying configure.sh


** Instructions **

1. Create a global project directory in your scratch (e.g., /scratch/Your_organims

2. Inside the global project directory, create a directory called "fastqs"
and a directory called "project.1".  Each run of the pipeline should have
an additional project directory ("project.2","project.3", etc), but
"fastqs" directory remains the same.

3. manually create a table called "fastqs_to_align.txt"
and save it in "project.1". [ Each time you rerun the pipeline you create such
a table and save it in the current project directory. ]
 
The  table is a tab-delimited table with 3 columns. Column 1 is the
names of the samples to be aligned and analyzed. These sample names
should be relatively short and should contain alpha-numeric characters
only, no special characters (exception, underscores is ok), and no 
whitespaces. Columns 2 and 3 are the names (not paths) of the fastq
files.

4. Create subdirectories inside "fastqs" directory with the
sample names in column 1 of "fastqs_to_align.txt". Copy the fastq files whose
names are in columns 2 and 3 of "fastqs_to_align.txt" to each of the 
appropriate subdirectories

5. Create a reference genome directory (anywhere, but possibly in your home)
and manually download genome reference.  You then need to index the files. 
You will need to run bwa index, samtools index, and GATK 
CreateSequenceDictionary.jar.

see: http://gatkforums.broadinstitute.org/discussion/1601/
how-can-i-prepare-a-fasta-file-to-use-as-reference

6. Copy the pipeline files (not the pipeline directory) to "project.1"

7. Modify configure.sh variables (additional information for each 
variable is found in the configure.sh itself) For FQTABLE, FQTOPDIR 
and REF variables, specify the absolute paths (from root) to "fastqs_to_align.txt" 
(e.g., /Path/To/fastqs_to_align.txt), the path to the "fastqs" directory 
(e.g., /Path/To/fastqs; note: do not add a trailing "/") and the path 
to your reference genome fasta file (e.g., /Path/To/reference.fasta). 
Also set the SAMPLEWISE_T_OPT and CHRWISE_T_OPT following the 
instructions specified in configure.sh. 

8. The GENOMEVERSION variable is necessary for snpEff.  It should match 
the (1) basename of your reference fasta (e.g., basename.fa or basename.fas). 
[[ in section 9, you will also see that GENOMEVERSION
should match the genome name in 2 different sections of the snpEff.config
file ]]


9. Please verify that your reference genome files use the same idenitifiers
(including exact same case) in the gff and fasta files. For example, if in
the gff chromosomes are referred to as, for example, "Chr1" they must also
be ">Chr1" in the fasta.


10. Execute the pipeline by cd'ing into "project.1" and type "bash master.sh".
This will generate the "alignments.1" directory where all outputs from
master.sh can be found.  You can monitor the progress using qstat -u youruserid.

Note: master.sh is a convenience method. It is possible to run each step consecutively
#using the qsub commands in master.sh. If do so, you can remove the -W options
#(while still honoring the ordre of scripts, e.g., snps cant be called until
#aln_pipeline is complete), but you will need to specify the -t options where applicable

11. After the pipeline has completed you must check the log files for errors. This is important
as hpc crashes frequently and for reasons out of your control (e.g., someone else 
crashes the node your job is running on).  

To check for errors, do the following:
In alignments.1, snps.1, snps.2, indels.1, indels.2 directories:
grep _ESTATUS_ *log | less  (check that all exit status of executed processes are zero)

or for only non-zero exit statuses,
grep _ESTATUS_ *log | grep -v -P "\s+0$" | less

grep _ERROR_ *log | less (check that no errors were found by the pipeline itself)
grep _END_ *log | less (check that each shell script completed)

*note: you should grep for errors separately in each of the mentioned directories.

Finally, if you would like to doublecheck the read counts in various versions of alignments
during aln_pipeline.sh, you can grep _COUNTS_ for each log file separately.  This is 
a way to track what is happening during the  removal of duplicates and the mates of duplicates
that are unmapped (and not marked as duplicates).

12. To re-run phase 1 of the pipeline (say, with different parameters), simply create "project.2"
inside your global project directory and repeat the above steps. Again, "fastqs" can
remain the same and if the fastq files are not different then the same "fastqs_to_align.txt"
file can be copied into "project.2"

13. To run the second phase of the analysis. Copy INSERTION_DETECTION.sh and INSERTION_DETECTION.pbs
into alignmnet.1 directory and run it by using following command: 

qsub -t INSERTION_DETECTION.pbs # t is the number of samples -1 in order to perform parallel proces. 

14. This part of the pipeline will perform quality control for mapping, coverage filtering 
and find where all the insertions are in the relationship to genes in the genome. 

** Notes **

(1) we are not trimming reads before alignment. Trimming has been shown to increase
the false positive rate

#General help on running PBS jobs:
#http://lnf.umich.edu/nnin-at-michigan/index.php/computation/running-jobs/

#General help on PBS dependency types (before, after, afterok, afternotok, etc.).  For more examples, see the qsub manpage or
#http://www.clusterresources.com/torquedocs21/commands/qsub.shtml

#The full range of dependency operaters is defined here:
#http://www.clusterresources.com/torquedocs21/commands/qsub.shtml

#Additional tips on dependency pbs scripting here:
#http://beige.ucs.indiana.edu/I590/node45.html

#**can set the slot limit on an array job using -t 0-99%25 (max 25 at a time in this example)
#for really big jobs

