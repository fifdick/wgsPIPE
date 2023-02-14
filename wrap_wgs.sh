#!/bin/sh 
#SBATCH --account=p94
#SBATCH --time=30:00:00
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

sID="censored"
key="01_blood"

input1="/cluster/projects/p94/fdi/fastqs/"$sID"_R1.fastq.gz"
input2="/cluster/projects/p94/fdi/fastqs/"$sID"_05_R2.fastq.gz"



mkdir -p /cluster/projects/p94/fdi/"$key"

outputpath="/cluster/projects/p94/fdi/"$key""
p="/cluster/projects/p94/fdi/scripts" #where wgs_pipe.sh script is located

chmod +x "$p"/wgs_pipe.sh
"$p"/wgs_pipe.sh -h > "$outputpath"/help.txt

echo "trimming"
date
if "$p"/wgs_pipe.sh -o $outputpath -i $input1 $input2 -sID $sID -t tr > "$outputpath"/trimming.txt ; then
echo "trimming successfull"
date
else 
echo "trimming failed"
exit 1
fi

date
echo "aligning"
if "$p"/wgs_pipe.sh -o $outputpath -sID $sID -c 16 -t aln > "$outputpath"/aligning.txt ; then
echo "aligning successfull"
date
else
echo "Aligning failed"
 exit 1
fi

date
echo "recalibrating"
if "$p"/wgs_pipe.sh -o $outputpath -sID $sID -i "$outputpath"/dedup_"$sID".bam -c 8 -t brec > "$outputpath"/recalibration.txt ; then 
echo "recalibrating successfull"
date
else 
echo "Recalibrating basescores  failed"
exit 1
fi

date
echo "variant calling"
if "$p"/wgs_pipe.sh -o   $outputpath -sID $sID -i "$outputpath"/recal_reads_"$sID".bam -c 8 -t varcall "$outputpath"/variantcalling.txt ; then 
echo "varcalling exit code 0"
date
else 
echo " variant calling failed "
exit 1
fi

exit 

