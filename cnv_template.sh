#!/bin/sh 
#SBATCH --account=p94
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --nodes=1 --exclusive
#SABTCH --ntasks-per-node=16

source /cluster/bin/jobsetup
set -o errexit

#setup tool dependencies
source /cluster/software/VERSIONS/CNVnator/setup.source
#reference data
ref="/cluster/shared/bioinformatics/reference-data/b37/GATK2.8/human_g1k_v37_decoy.fasta"

##########INPUT FROM JS WILL COME AFTER THIs



mkdir -p "$outdir"

#CNVnator steps:
#first step
echo "starting CNV pipeline with CNVnator on sample:$sample_name"
echo "results will be written to $outdir"
echo "CNVnator 1st step:"
if cnvnator -genome GRCh37 -root "$outdir"/out.root -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y -tree "$input_bam" -unique; then 
echo "Moving on to 2nd step"
else "first step of CNVnator failed"
exit 1;
fi
if cnvnator -genome GRCh37 -root "$outdir"/out.root -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y -his 100 -d /cluster/projects/p94/fdi ; then
echo "Moving on to 3rd step"
else "second step of CNVnator failed"
exit 1;
fi
if cnvnator -root "$outdir"/out.root -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y -stat 100 ; then
echo "Moving on to 4th step"
else "third step of CNVnator failed"
exit 1;
fi
if cnvnator -root "$outdir"/out.root -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y -partition 100 ; then
echo "Moving on to last step of CNVnator"
else echo " fourth step of CNVnator failed"
exit 1
fi
if cnvnator -root "$outdir"/out.root -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y -call 100 > ""$outdir"/CNVnator/"$sample_name".txt" ; then 
echo "sucessfully caled CNVs for sample "$sample_name""
else
echo "failed in last step of CNVnator"
exit 1
fi

#ERDS
if /cluster/software/VERSIONS/CNVnator/erds/erds1.1/erds_pipeline.pl -b "$input_bam" -v "$input_vcf" -o "$outdir/ERDS" -r "$ref" -sd b37 -name "$sample_name" ; then
echo "erds performed sucessfully"
exit
else 
echo "erds failed on "$sample_name""
exit 1
fi

