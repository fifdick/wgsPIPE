#!/bin/sh 
#SBATCH --account=p94
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --nodes=1 --exclusive
#SABTCH --ntasks-per-node=16

source /cluster/bin/jobsetup
set -o errexit


#TODO
#create variables for ids, loop through samples
#variables for tool versions and dependecies
#set versions with file??!
#check files, versions -tools etc.
#blabla


mkdir -p /cluster/projects/p94/fdi/gvcf_results/8_sample_run/

outputpath="/cluster/projects/p94/fdi/gvcf_results/8_sample_run"

p="/cluster/projects/p94/fdi/scripts" #where wgs_pipe.sh script is located

chmod +x "$p"/wgs_pipe.sh

vcfs=":/cluster/projects/p94/fdi/03_brain/69_15.g.vcf:/cluster/projects/p94/fdi/03_blood/10_08.g.vcf:/cluster/projects/p94/fdi/04_brain/284_10.g.vcf:/cluster/projects/p94/fdi/04_blood/159_09.g.vcf:/cluster/projects/p94/fdi/gvcf_results/jointcalls/multisample.g.vcf"

now=$(date)

#date
#echo "merging vcfs with gatk CombineVariants and joint calling on merged file afterwards"
#if "$p"/wgs_pipe.sh -o $outputpath -i "$vcfs" -c 4 -t joint_call > "$outputpath"/merge_joint_log_"$now".txt ; then 
#echo "merging exit code 0"
#date
#else 
#echo "merging failed "
#exit 1
#fi

#now=$(date)
#if "$p"/wgs_pipe.sh -o $outputpath -i "$outputpath/multisample.g.vcf" -c 4 -t joint_call > "$outputpath"/mergelog_"$now".txt ; then
#echo "succesfully performed joint_calling"
#exit
#else
#echo "joint_calling failed"
#exit  1
#fi
 
#now=$(date)
#if "$p"/wgs_pipe.sh -o "$outputpath" -i ""$outputpath"/joint_call.vcf" -c 4 -VQSR_tr 90.97 -t var_recal > "$outputpath"/var_recal_log_"$now".txt ; then
#echo "succesfully performed variant recalibration"
#else
#echo "variant_recailbration failed"
#exit  1
#fi
    

now=$(date)
if "$p"/wgs_pipe.sh -o "$outputpath" -i ""$outputpath"/recalibrated_variants.vcf" "08_14:189_05" -c 8 -t filter_var > "$outputpath"/filter_log_"$now".txt ; then
echo "succesfully performed variant filtering 1"
else
echo "variant filtering 1 failed"
exit  1
fi


if "$p"/wgs_pipe.sh -o "$outputpath" -i ""$outputpath"/recalibrated_variants.vcf" "29_15:25_05" -c 8 -t filter_var >> "$outputpath"/filter_log_"$now".txt ; then
echo "succesfully performed variant filtering 2"
exit
else
echo "variant filtering 2  failed"
exit  1
fi

