
#!/bin/bash

#General help on running PBS jobs:
#http://lnf.umich.edu/nnin-at-michigan/index.php/computation/running-jobs/

#General help on PBS dependency types (before, after, afterok, afternotok, etc.).  For more examples, see the qsub manpage or
#http://www.clusterresources.com/torquedocs21/commands/qsub.shtml

#The full range of dependency operaters is defined here:
#http://www.clusterresources.com/torquedocs21/commands/qsub.shtml

#Additional tips on dependency pbs scripting here:
#http://beige.ucs.indiana.edu/I590/node45.html

#**can set the slot limit on an array job using -t 0-99%5 (max 5 at a time)

#Instructions:
#1. manually create "fastqs_to_align.txt"
#2. manually download genome reference and run indexing script
#3. configure pipeline to run by (1) adding # to skips to skip in this script
#                                (2) configure the configure.sh script 

#execute configure.sh to define environment
if [ -f configure.sh ]; then
        source configure.sh
else
        echo _ERROR_ configure.sh not found
        exit 1
fi

cd alignments.1

#Sample-level alignment
ALN=`qsub -t ${SAMPLEWISE_T_OPT} aln_pipeline.pbs`

