#!/usr/bin/perl

use strict;
use warnings;

use 5.010;

#cnv_pipe_template="CT"
my $CT;
#job_script_cnv_sample_specific="JS"
my $JS;

my $vcf_dir="/cluster/projects/p94/fdi/gvcf_results/8_sample_run";

my $CNV_out="/cluster/projects/p94/fdi/CNV";
system("mkdir -p $CNV_out");

foreach my $file (@ARGV) {
if(-e $file )
{
#the file that will be th einput to CNVnator (in the jobs script that needs to be generated
#We thus need the file itself AND the name to add the string to the jobscriptname

my $recal_read=$file;
say "Starting to process the following file:";
say "$recal_read";

#extracting the samplename from the file path (only possible assuming the it lies in this level of dir hirachy:
#clutser > projects > p94 > fdi > recals (or some other foldername) > file
my @split_path = split /[\/.]/, $recal_read;
my $sample_name = $split_path[6];

my $vcf_input;
opendir(DIR, $vcf_dir) or die "could not open $vcf_dir";
my @vcf_files= grep { /.*${sample_name}{1}.*\.vcf$/ } readdir(DIR);
closedir(DIR) or die "died trying to close $vcf_dir";

if ( scalar @vcf_files != 1 )
{
say "Found more than one match for vcf filename in ${vcf_dir}, probably something wrong with filename (should be only 2 ids of the two samples represented in the vcf)\n or sth wen wrong during grep regex" ;
}
else 

{
$vcf_input="$vcf_dir/$vcf_files[0]";
say "using $vcf_input as input for ERDS pipeline";
}


say "Now generating a job_script for $sample_name";

#create a job script with the sample name in it's filename string
open(JS,">/cluster/projects/p94/fdi/scripts/jobscripts/$sample_name.sh") or die "ERROR:failed to create job script for sample: $sample_name";
open(CT,'/cluster/projects/p94/fdi/scripts/cnv_template.sh') or die "ERROR:CNV_template file not found";

while(<CT>)
{
	if ( $_ =~ /^#{10}/)
	{	
		print JS 'sample_name="' . $sample_name . "\"\n";
		print JS 'input_bam="' . $recal_read . "\"\n";
		print JS 'outdir="' . "$CNV_out/$sample_name" . "\"\n";
		print JS 'input_vcf="' . "$vcf_input" . "\"\n";
		next;
	} 
	else 
	{print JS $_;}
}



close JS or die "cant close file jobscript";
close CT or die "cant close cnvpipe template";

system("sbatch /cluster/projects/p94/fdi/scripts/jobscripts/$sample_name.sh")
}} 
