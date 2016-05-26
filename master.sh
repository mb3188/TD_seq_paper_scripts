
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
#4. modify pindel.bamcfg
#5. modify pindel.chr2run (names of chromomosomes must be same as reference genome fasta identifiers)

#execute configure.sh to define environment
if [ -f configure.sh ]; then
        source configure.sh
else
        echo _ERROR_ configure.sh not found
        exit 1
fi

mkdir alignments.1
cp aln_pipeline.sh \
aln_pipeline.pbs \
merge_bams.sh \
merge_bams.pbs \
coverage_depth.sh \
coverage_depth.pbs \
indel_realignment.sh \
indel_realignment.pbs \
callable.sh \
callable.pbs \
callable_merge.sh \
callable_merge.pbs ./alignments.1

mkdir alignments.1/snps.1
cp gatk_pipeline.sh \
gatk_pipeline.pbs \
snpEff_pipeline.sh \
snpEff_pipeline.pbs \
snpqc.pbs \
snpqc.sh \
snpqc.R \
variantrecalibrator.sh \
variantrecalibrator.pbs \
./alignments.1/snps.1

mkdir alignments.1/snps.2
cp mpileup.pbs \
mpileup.sh \
snpEff_pipeline.sh \
snpEff_pipeline.pbs \
snpqc.pbs \
snpqc.sh \
snpqc.R ./alignments.1/snps.2

mkdir alignments.1/indels.1
cp gatk_indels.sh \
gatk_indels.pbs ./alignments.1/indels.1

mkdir alignments.1/sv.1
cp pindel.sh \
pindel.pbs \
pindel.bamcfg \
pindel.chr2run \
pindel2vcf.sh \
pindel2vcf.pbs ./alignments.1/sv.1

mkdir alignments.1/sv.2
cp delly.sh \
delly.pbs alignments.1/sv.2

mkdir alignments.1/snps.1/filtered.1
cp filter_vcf.pbs \
filter_vcf.sh ./alignments.1/snps.1/filtered.1

mkdir alignments.1/snps.2/filtered.1
cp filter_vcf.pbs \
filter_vcf.sh ./alignments.1/snps.2/filtered.1

cd alignments.1

#Sample-level alignment, fix-mates, re-alignment, de-dupping [ < 24 hrs chlamy, ~24 hrs dates ]
ALN=`qsub -t ${SAMPLEWISE_T_OPT} aln_pipeline.pbs`

# note that I put %30 here in order to not allow ot do alignment onmore than 30 genomes at the time so the cluster does not go crazy 

#delly detects deletions, but does not use merged bams
cd sv.2
DELLY=`qsub -W depend=afterokarray:${ALN} -v SBAM=bwa_msdr.bam delly.pbs`
cd ../

CALLABLE=`qsub -W depend=afterokarray:${ALN} -t ${SAMPLEWISE_T_OPT}%30 -v SBAM=bwa_msdr.bam callable.pbs`
CALLABLEMERGE=`qsub -W depend=afterokarray:${CALLABLE} -v SBAM=bwa_msdr.bam callable_merge.pbs`

# note that I put %30 here in order to not allow ot do alignment onmore than 30 genomes at the time so the cluster does not go crazy  

#Merge bams 
MERGE=`qsub -W depend=afterokarray:${ALN} -v SBAM=bwa_msdr.bam merge_bams.pbs`

#get coverage [ 35 hrs on merged chlamy genomes ]
COVERAGE=`qsub -W depend=afterok:${MERGE} -v MBAM=bwa_msdr_M.bam coverage_depth.pbs`


######################

#  I NEED TO DO GLOBAL REALIGNMENT AND I AM DEFINING THIS VARIABLE HERE

GLOBALREALIGNMENT=1



#Global indel realignment [ 40+ hrs on chlamy genomes ]
if [ $GLOBALREALIGNMENT -eq 1 ]; then

	REALN=`qsub -W depend=afterok:${MERGE} -v MBAM=bwa_msdr_M.bam indel_realignment.pbs`

        cd snps.1
        CALLSNPS=`qsub -W depend=afterok:${REALN} -v MBAM=bwa_msdr_MR.bam gatk_pipeline.pbs`
	SNPQC=`qsub -W depend=afterok:${CALLSNPS} -v INVCF=bwa_msdr_MR.vcf snpqc.pbs`
	SNPEFF=`qsub -W depend=afterok:${CALLSNPS} -v INVCF=bwa_msdr_MR.vcf snpEff_pipeline.pbs`

	cd ../snps.2
	MPILEUP=`qsub -W depend=afterok:${REALN} -v MBAM=bwa_msdr_MR.bam mpileup.pbs`
	MPILEUPSNPQC=`qsub -W depend=afterok:${MPILEUP} -v INVCF=bwa_msdr_MR.vcf snpqc.pbs`

	cd ../sv.1
	PINDEL=`qsub -W depend=afterok:${REALN} -t ${CHRWISE_T_OPT} pindel.pbs`
	PINDEL2VCF=`qsub -W depend=afterokarray:${PINDEL} pindel2vcf.pbs`	

	cd ../indels.1
	CALLINDELS=`qsub -W depend=afterok:${REALN} -v MBAM=bwa_msdr_MR.bam gatk_indels.pbs`
	
	
elif [ $GLOBALREALIGNMENT -eq 0 ]; then

        cd snps.1
        CALLSNPS=`qsub -W depend=afterok:${MERGE} -v MBAM=bwa_msdr_M.bam gatk_pipeline.pbs`
        SNPQC=`qsub -W depend=afterok:${CALLSNPS} -v INVCF=bwa_msdr_M.vcf snpqc.pbs`
	SNPEFF=`qsub -W depend=afterok:${CALLSNPS} -v INVCF=bwa_msdr_M.vcf snpEff_pipeline.pbs`

	cd ../snps.2
        MPILEUP=`qsub -W depend=afterok:${MERGE} -v MBAM=bwa_msdr_M.bam mpileup.pbs`
	MPILEUPSNPQC=`qsub -W depend=afterok:${MPILEUP} -v INVCF=bwa_msdr_M.vcf snpqc.pbs`

	cd ../sv.1
	PINDEL=`qsub -W depend=afterok:${MERGE} -t ${CHRWISE_T_OPT} pindel.pbs`
	PINDEL2VCF=`qsub -W depend=afterokarray:${PINDEL} pindel2vcf.pbs`

	cd ../indels.1
	CALLINDELS=`qsub -W depend=afterok:${MERGE} -v MBAM=bwa_msdr_M.bam gatk_indels.pbs`
else
	echo _ERROR_ GLOBALREALIGNMENT boolean invalid in master.sh
	exit 1
fi


MYID=`whoami`
qstat -u ${MYID} > master.log
