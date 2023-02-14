#!/bin/sh
#SBATCH --account=p94
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --nodes=1 --exclusive
#SABTCH --ntasks-per-node=16

source /cluster/bin/jobsetup
set -o errexit

cmd="$1"
param1="$2"
param2="$3"


if [ -z $cmd ] ; then 
echo "no command or script specified"
exit 1;
else
	if $cmd $param1 $param2 ; then
	echo "command exit with code 0"
	else 
	echo "command failed"
	exit 1
	fi
exit 
fi
