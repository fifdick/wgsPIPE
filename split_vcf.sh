#!/bin/sh

module load bcftools

vcf_dir="$1"
out_dir="$2"

if [ -z "$out_dir" ] ; then
outdir="$vcf_dir"
fi

echo "splitting files from dir : "$vcf_dir" to "$out_dir""
for file in "$vcf_dir"/*.vcf*; do
	echo "now splitting file "$file""
	for sample in `bcftools query -l $file`; do
		echo "looking at sample  "$sample""
		bcftools view -c1 -Ov -s $sample -o "$out_dir"/"$sample".vcf $file
	done
done
echo "done"
