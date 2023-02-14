#!/bin/sh
# exit on error
set -o errexit 
# force variables to be set
#set -o nounset #CANNOT USE THIS YET BECAUSE FOR ALN THE INPUT CAN BE UNSET AND WILL BE SET LATER...
# force wildcards to throw error
shopt -s failglob

#TODO INCLUDE VARIABLES FOR VERSIONS....
#TODO CHECK PATHS AND TOOLS 

#check PATHS
### Tool Paths ##
#define constant variables that wont be changed throuhout the script
#readonly java="/cluster/software/VERSIONS/java/jdk1.8.0_112/bin/java"
#readonly trimmomatic="/cluster/shared/bioinformatics/src/Trimmomatic-0.33/trimmomatic-0.33.jar"
#readonly picard="/cluster/shared/bioinformatics/src/picard-tools-1.129/picard.jar"
#readonly gatk="/cluster/shared/bioinformatics/src/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar"

### FILE PATHS ###

#readonly adapters="/cluster/shared/bioinformatics/src/Trimmomatic-0.33/adapters/TruSeq3-PE-2.fa"
#readonly reference="/cluster/shared/bioinformatics/reference-data/b37/GATK2.8/human_g1k_v37_decoy.fasta"
#readonly known_snps="/cluster/shared/bioinformatics/reference-data/b37/GATK2.8/dbsnp_138.b37.vcf"
#readonly known_indels="/cluster/shared/bioinformatics/reference-data/b37/GATK2.8/Mills_and_1000G_gold_standard.indels.b37.vcf"


#Help message printed if called with > wgs_pipe.sh -h 
usage() {

echo "Preprocessing pipeline for Whole Genome Sequencing Data (PE)"
echo "" 
echo "Usage: $0 -o <OUTPUTDIR> -i <input_1> <input_2> -c <#threads> -sID <sampleID> -t <\"tr\"/\"aln\"/\"brec\"/\"varcall\"/\"varcall_compare\" > "
echo ""
echo "!!!!!!!!!!!!!!!Please keep the order in the right way, as defined above !!!!!!!!!!!!!!!!!!!!!!"
echo ""
echo "Mandatory arguments:"
echo " -o, --output <path>    Specifies output Path."
echo " -i, --input  <input_1> <Input_2> "
echo "		For trimming: specifiy forward and reverse reads in 2 input files."
echo "		For alignment , no files need to be provided if trimming has been done before alignment and trimmed reads are stored in -o <OUTPUTDIR> (If trimmed reads for alignment are not provided, ResultDirs will be searched for trimmed reads with the file extension \"*_1P.fq.gz\" and \"*_2P.fq.gz\")"
echo "		For base recalibration only one Input file needed (trimmed, aligned, sorted, indexed and deduped)"
echo "		For variant calling only one Input file needed (recalibrated reads)"
echo " 		For variant calling comparing alleles, two Input files are required: Inputfile1 should be the bam file and Inputfile 2 should be a vcf file for comparing alleles."
echo " 		For joint calling an prior merging (if multiple g.vcfs) please specify gvcfs as a string delimited by colon and a leading colon. Like so"
echo "			:/path/first.g.vcf:/path/second.g.vcf:/path/andsoon.gvcf	"
echo " 		For joint calling on a already merged multisample.g.vcf just specify the file as a normal input file"
echo "		For variant recalibration specify one inputfile (.vcf) and optionally VQSLOD threshold wit the option -VQSR_tr otherwise 99.0 will be used"
echo "-t, --tool <\"tr\"/\"aln\"/\"brec\"/\"varcall\"/\"varcall_compare\"/\"joint_call\"/\"var_recal\" > Specifies tool."
echo "		Possible arguments: -tr for trimming , -aln for alignment , -brec for base recalibration "
echo "				    -varcall for variant calling in discovery mode"
echo " 				    -varcall_compare for variant calling only those alleles provided in extra inputfile" 
echo "   			    -joint_call for joint calling of one or multiple g.vcf files"
echo "				    -var_recal for variant recalibration according to VQSR gatk's best practice"
echo " " 
echo "OPTIONS"
echo " -VQSR_tr , --VQSR_threshold setting the VQSLOD threshold that will be used during variant recalibration (default = 99.0)"
echo " -sID , --sampleID Id info of sample will be incorperated to ReadGroup info druing alignment"
echo "			 "pseudo_id" will be used if no id provided."
echo " -c , --cpus Number of threads used for each command. Default is 4"
#echo " *  If neither \"-tr\" option, nor \"-aln\" option -> both will be executed consecutively"
}

#############  PARSE CONFIG FILE   #################################
###############################################################
set_dependencies() {
echo " The following tools/software are required: "
echo " JAVA, trimmomatics, picard, gatk"

tools=($(echo 'cat //tool/*/@name' | xmllint --shell config_wgs.xml | awk -F\" 'NR % 2 == 0 {print $2; }'))
versions=($(echo 'cat //tool/*/@version' | xmllint --shell config_wgs.xml | awk -F\" 'NR % 2 == 0 {print $2; }'))
locations=($(echo 'cat //tool/*/@location' | xmllint --shell config_wgs.xml | awk -F\" 'NR % 2 == 0 {print $2; }'))


echo "You have provided the following paths in your config file (config_wgs.xml)"
for ((i=0;i<${#tools[@]};++i)); do
	software=${tools[i]}
	location=${locations[i]}
	version=${versions[i]}
	path=${location}${version}	
	eval $software=\$path
	echo " ${software}=${path} "
		
done

files=($(echo 'cat //file/*/@name' | xmllint --shell config_wgs.xml | awk -F\" 'NR % 2 == 0 {print $2; }'))
locations=($(echo 'cat //file/*/@location' | xmllint --shell config_wgs.xml | awk -F\" 'NR % 2 == 0 {print $2; }'))
echo " You have specified the following file paths in you config file: "
for ((i=0;i<${#files[@]};++i)); do
        file_name=${files[i]}
        location=${locations[i]}
        eval $file_name=\$location
        echo " ${file_name}=${location} "

done



}
#Test commandline parsing
parsing() {
if [ -z $1 ] || [ -z $2 ]; then

	echo "One or both inputfile paths are undefined"
        exit 1
elif [ ! -f "$1" ]; then
	echo "Inputfile $1 doesn't exist"
	exit 1
elif [ ! -f "$2" ]; then
	echo "Inputfile $2 doesn*t exist"
	exit 1
fi 
echo "You entered the following input files:${1} and ${2}"
}


load_modules() {
### LOAD MODULES ###
module purge
arr=("R" "samtools" "bwa" "gatk")

R_version="/$(echo 'cat //R/*/@version' | xmllint --shell config_wgs.xml | awk -F\" 'NR % 2 == 0 {print $2 }')"
samtools_version="/$(echo 'cat //samtools/*/@version' | xmllint --shell config_wgs.xml | awk -F\" 'NR % 2 == 0 {print $2 }')"
bwa_version="/$(echo 'cat //bwa/*/@version' | xmllint --shell config_wgs.xml | awk -F\" 'NR % 2 == 0 {print $2 }')"


arr_v=("$R_version" "$samtools_version" "$bwa_version")
for (( i=0; i < ${#arr_v[@]}; i++)); do
	if [ -z ${arr_v[$i]} ]; then 
		echo "you did not specify a version for one of the modules in your config_wgs.xml!--> loading default"
	fi
done

for (( i=0; i < ${#arr[@]}; i++)); do
	echo "... loading module ${arr[$i]} ... "
	echo "module load ${arr[$i]}${arr_v[$i]}"
	module load ${arr[$i]}${arr_v[$i]}
		if [ $? -ne 0 ]; then
                        echo "loading module ${arr[$i]} failed"
                        exit 1
                fi	
	echo "from this path:"
	#get path from "module show <tool>" output by redirecting output to stdout (before stderr)
	module show ${arr[$i]} 2>&1 | grep -m1 PATH | awk -F '\t' '{print $2}'
done
}

set_vars () {

#DEFAULT VARIABLES
if [ -z "$VQSR_tr" ]; then
VQSR_tr="99.0"
fi

sample_name="pseudo_samplename"                                                                                   
lib_id="pseudo_libname"
platform="ILLUMINA" 
#check if sample id was provided by user and has been parsed 
# if not, use "pseudo_id"
if [ -z "$sID" ]; then 
	echo "You did not provide a sample ID"
	sID="pseudo_id" #could be set to different value  by user with commandline arg -sID
	readgroups="@RG\tID:"$sID"\tSM:"$sample_name"\tLB:"$lib_id"\tPL:"$platform"\tPU:"$sID""
else 
	readgroups="@RG\tID:"$sID"\tSM:"$sID"\tLB:"$lib_id"\tPL:"$platform"\tPU:"$sID""
fi
#If # of cpus was not specified by user, set it to 4 by default
if [ -z "$cpus" ]; then
	cpus="8" #could be set to different value  by user with commandline arg -c
fi

}

#function to check if variables have been set to prevent usage of unset variables in functions later on 
check_vars () {

varlist="$1"
command_name="$2"

for i in "${varlist[@]}"
do 
	if [ -z "$i" ]; then
	echo "one of the variables: "$i" used in the "$command_name" command are not set correctly. Please check the script (function \" "$command_name" () \" ) "
	exit 1
	fi
done
	

}

#################################################
### TRIMMING ###
#################################################

trimming() {

echo "... Trimming reads"
#Set input variables which have been parsed from commandline and provided to function during function call 

fastqInput_1=$1
fastqInput_2=$2

#CHECK IF ALL VARIABLES HAVE BEEN SET CORRECTLY
declare -a vars=("$java" "$trimmomatic" "$resultDir" "$fastqInput_1" "$fastqInput_2" "$adapters" "$resultDir")
check_vars $vars "trimming" #function will abort script if one of those variables are unset

#start trimming
date
if $java -jar $trimmomatic PE -phred33 -baseout "$resultDir"/"$sID".fq.gz "$fastqInput_1" "$fastqInput_2" ILLUMINACLIP:$adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36; then
	echo "trimming was successfull"
	date
else 
	echo "Non-zero exit status. Trimming failed"
	exit 1
fi

}
 
##################################################
###ALIGNMENT###
##################################################

alignment() {

#Set input variables which have been parsed from commandline or were found in the specified output directory
#TODO cannot yet handle case when multiple files  found in dir

alnInput_1=$1
alnInput_2=$2

#Set output variable
sorted_bam=""$resultDir"/sorted_"$sID".bam"
#CHECK IF ALL VARIABLES HAVE BEEN SET CORRECTLY
declare -a vars=("$cpus" "$reference" "$alnInput_1" "$alnInput_2" "$sorted_bam" "$resultDir" "$picard" "$java" "$readgroups")
check_vars $vars "alignment" #function will abort script if one of those variables are unset

#start alignment process
echo "... aligning reads (PE) with BWA mem, converting to bam and sorting with samtools ...."
echo "The sampleID you entered  will be incorperated into the readgroup info like so:"
echo "$readgroups"
date
#set pipefail so that even if last command in pipe exited with 0, pipe with exit with 1 if error occured inbetween
#set -euxo pipefail DOES NOT WORK HOW IT SHOULD:::THROWS ERROR OF UNSET VAR BEFORE VAR GETS SET IN SEQUENCE OF PIPE?!??!?!?!??!?
#TODO test if pipefail setting works correctly?


#check if referece file exits
if [ -f "$reference" ]; then # "-R" in the bwa mem command for adding readgroups!
	bwa mem -t "$cpus" -R "$readgroups" "$reference" "$alnInput_1" "$alnInput_2" | samtools view -b -@ "$cpus" - | samtools sort -O BAM -o "$sorted_bam" -@ "$cpus" -
	if [ "$?" -ne 0 ]; then echo "bwa command failed"; exit ;
	else
		 echo "aln sucessfull"
		 date
	fi
else
	echo "reference not found"
	exit 1
fi

#output file written to same location as "$sorted_bam"
echo "... indexing aligned bam file..."
if [ -f "$sorted_bam" ]; then
	samtools index $sorted_bam
	if [ "$?" -ne 0 ]; then echo "Failure during indexing of the sorted bam file"; exit 1; fi
else echo "no sorted bam file found to be indexed."
fi

### REMOVING DUPLICATES ###
echo "... removing duplicates with Picard..."
date
if [ -f "$sorted_bam" ]; then
	$java -Xmx4g -Djava.io.tmpdir="$resultDir" -jar $picard MarkDuplicates I="$sorted_bam" O=""$resultDir"/dedup_"$sID.bam"" CREATE_INDEX=true REMOVE_DUPLICATES=TRUE M=""$resultDir"/output.metrics"
	if [ "$?" -ne 0 ]; then echo "Removing/Marking Duplicates failed"; exit 1;
	else
	echo "dedup sucessfull"
	date
	fi 
else "no bam file found to remove/mark duplicates"
fi
}

##########################################################
### BASE RECALIBRATION ###
##########################################################

base_recal () {

#Base recalibration with gatk

#Set Variables
deduped_bam="$1" #provided input file specified by user and submitted through funciton call during commandline parsing
br_out_1=""$resultDir"/recal_data_"$sID".table" 
br_out_2=""$resultDir"/post_recal_data.table" #NOT USED BECAUSE GGPLOT NOT WORKING ON TSD
plot=""$resultDir"/recal_plot.pdf" #NOT USED BECAUSE GGPLOT NOT WORKING ON TSD YET
recal_out=""$resultDir"/recal_reads_"$sID.bam""

#CHECK IF ALL VARIABLES HAVE BEEN SET CORRECTLY
declare -a vars=("$java" "$known_snps" "$known_indels" "$reference" "$deduped_bam" "$br_out_1" "$br_out_2" "$gatk", "$plot", "$resultDir")
check_vars $vars "base_recal" #function will abort script if one of those variables are unset

#start base recalibration process
echo " Step one of base recalibration. Calculating covariation data. Saving $br_out_1"
date

if [ -f $known_snps ] && [ -f $known_indels ] && [ -f $reference ]; then #make sure files exist and the path provided is correct
	$java -jar $gatk -T BaseRecalibrator -nct "$cpus" -R "$reference" -I "$deduped_bam" -knownSites "$known_snps" -knownSites "$known_indels" -o "$br_out_1"
	if [ "$?" -ne 0 ]; then 
		echo "first step of base recalibration (covariant calculations) failed"
		exit 1
#	else #NOT EXECUTED BECAUSE GGPLOT LIBRARY FOR R IS NOT AVAILABLE ON TSD
#		echo " Second step and third of Base Recalibration: calc covariant after calibration for plotting (to compare later on) "
#		#$java -jar $gatk -T BaseRecalibrator -R $reference -I $sortdidxd_bam -knownSites $known_snps -knownSites $known_indels -BQSR $br_out_1 -o $br_out_2
##		$java -jar $gatk -T AnalyzeCovariates -R $reference -l DEBUG -before $br_out_1 -after $br_out_2 -plots $plot
#		if [ "$?" -ne 0 ]; then
#			echo "Second step of Base recalibration failed"
#   		fi
	fi
	echo " Last step of Base Recalibration: applying calculations to input file"
	$java -jar $gatk -T PrintReads -nct "$cpus" -R $reference -I $deduped_bam -BQSR $br_out_1 -o $recal_out
	if [ "$?" -ne 0 ]; then
		echo " last step of base recalibration failed"
		exit 1
	fi
	echo "Base Realibration completed sucessfully"
	date
	exit
	
else
	echo "fetching $known_snps and $known_indels failed. Can't find files."
	exit 1
fi
}

############################################################
### VARIANT CALLING ####
############################################################

#TODO
varcall_discovery() {

#Calling varicants with gatk haplotype caller and default parameters
processed_bam="$1" #bam file provided by user
varcall_out=""$resultDir"/"$sID".g.vcf"

#check that no variables that will be used are unset
declare -a vars=("$java" "$gatk" "$processed_bam" "$reference" "$varcall_out" "$resultDir")
check_vars $vars "varcall" #function will abort script if one of those variables are unset

echo "Calling varinats with gatk's Haplotype Caller (.g.vcf)"
#Make sure reference file and input file exist
if [ -f "$reference" ] && [ -f "$processed_bam" ]; then
	if "$java" -jar "$gatk" -T HaplotypeCaller -nct "$cpus" -R "$reference" -I "$processed_bam" --genotyping_mode DISCOVERY --emitRefConfidence GVCF -stand_call_conf 30 -o "$varcall_out"; then
		echo "successfully called variants. Output can be found in "$varcall_out""
		date 
	else
		echo "variant calling failed"
		exit 1
	fi
else
	echo "the reference file: "$reference" or the input file: "$processed_bam" doesnt exit or couldnt be found" 		
	exit 1
fi	
}

varcall_given_alleles() {

#variant calling comparing to given vcf

processed_bam="$1" #bam file for which variants should be called
alleles="$2" # only variants that were provided by this file are being called
varcall_out=""$ResultDir"/"$sampleID".vcf"

#make sure nor variables are unset
declare -a vars=("$java" "$gatk" "$processed_bam" "$reference" "$varcall_out" "$alleles")
check_vars $vars "varcall" #function will abort script if one of those variables are unset

echo "Calling varinats with gatk's Haplotype Caller (only looking for variants provided in vcf file (genotyping_mode=GENOTYPE_GIVEN_ALLELE)"
if [ -f "$reference" ] && [ -f "$processed_bam" ] && [ -f "$alleles" ]; then
	if [ "$java" -jar "$gatk" -T HaplotypeCaller -nt "$cpus" -R "$reference" -I "$processed_bam" -alleles "$alleles" --GENOTYPE_GIVEN_ALLELES -stand_emit_conf 10 -stand_ecall_conf 30 -o "$varcall_out" ]; then
		echo "successfully called variants. Output can be found in "$varcall_out""
		date
	else
		echo "variant calling failed"
		exit 1
	fi
else
	echo "the reference file: "$reference", input file: "$processed_bam" or the file with given alleles "$alleles" doesnt exit or couldnt be found" 		
	exit 1
fi	
}
#############################################################
##########JOINT CALLING#####
#############################################################

joint_calling () {

file_string="$1"
mergeflag="$2"

#check if individual vcf files exist, exit else
#TODO: merged with id of inputs in filename string

merge_out=""$resultDir"/multisample.g.vcf"
joint_call_out=""$resultDir"/joint_call.vcf"

if [ "$mergeflag" -eq "1" ]; then 
	gatk_variant_string=$(sed 's/:/ --variant /g' <<< "$file_string")

	##make sure no variables are unset
	declare -a vars=("$resultDir" "$java" "$gatk" "$file_string" "$gatk_variant_string" "$reference" "$merge_out")
	check_vars $vars "jointcalling" #function will abort script if one of those variables are unset
	merge_command=""$java" -jar "$gatk" -T CombineGVCFs -R "$reference" -o "$merge_out""
	merge_command="$merge_command""$gatk_variant_string"
	echo "$merge_command"
	date
	echo "merging multiple g.vcfs to one"
	if [ -f "$reference" ] ; then
        	if eval "$merge_command" ; then
		echo "provided g.vcf files have been merged to "$merge_out""
		else
		echo "merging failed"
		fi
	else echo "reference file doesnt exist"
	fi
	date
	echo " joint genotyping"
        date
        if "$java" -jar "$gatk" -T GenotypeGVCFs -R "$reference" -o "$joint_call_out" -nt "$cpus" --useNewAFCalculator --variant "$merge_out" ; then
                echo "successfully performed joint genotyping on provided file "$merge_out", ouput can be found in "$joint_call_out""
                exit
        else
                echo "joint genotyping failed"
                exit 1
        fi
        date

else
		
	echo " joint genotyping using provided inputfile "$file_string" (skipped merging)"
	date
	if "$java" -jar "$gatk" -T GenotypeGVCFs -R "$reference" -o "$joint_call_out" -nt "$cpus" --useNewAFCalculator --variant "$file_string" ; then 
		echo "successfully performed joint genotyping on provided file "$merge_out", ouput can be found in "$joint_call_out""
		exit
	else
		echo "joint genotyping failed"
		exit 1
	fi
	date
fi
 
}

###########################################################
### Variant Recalibration VQSR ###
###########################################################
var_recal() {

multisample_in="$1"
tranche_threshold="$2"
vqsr_out_1=""$resultDir"/recal_snps_raw_indels.vcf"
vqsr_out_2=""$resultDir"/recalibrated_variants.vcf"

####STEP 1A: SNPS####
###################


snp_out_1=""$resultDir"/multisample.snp.model"
snp_out_2=""$resultDir"/multisample.snp.model.tranches"
snp_out_3=""$resultDir"/multisample.snp.model.plots.R"

declare -a vars=("$resultDir" "$java" "$gatk" "$multisample_in" "$hapmap" "$reference" "$hapmap" "$snps_1kg" "$known_snps" "$snp_out_1" "$snp_out_2" "$cpus" "$snp_out_3" "$vqsr_out_1")
check_vars $vars "variant recalibraiton I" #function will abort script if one of those variables are unset

if [ ! -f "$reference" ] ||  [ ! -f "$snps_1kg" ] || [ ! -f "$hapmap" ] || [ ! -f "$omni" ] || [ ! -f "$known_snps" ] ; then
echo "some of the specified reference file could not be found, please check paths in config file"
exit 1
fi

if
"$java" -jar "$gatk" -T VariantRecalibrator \
-R "$reference" \
--input "$multisample_in" \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 "$hapmap" \
-resource:omni,known=false,training=true,truth=true,prior=12.0 "$omni" \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 "$snps_1kg" \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "$known_snps" \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
--mode SNP \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
--recal_file "$snp_out_1" \
--tranches_file "$snp_out_2" \
-nt "$cpus" \
-rscriptFile "$snp_out_3" ; then
echo "first step of variant recalibration: calculation of snp recalibration, was successfull"
else
echo "first step of variant recalibration: calculation of snp recalibration, failed!"
exit 1
fi

####STEP 1B: SNPS####
###################
#TODO maybe use ts_filter_level 90.0 to exclude more FP

if
"$java" -jar "$gatk" -T ApplyRecalibration \
-R "$reference" --input "$multisample_in" --mode SNP --ts_filter_level "$tranche_threshold" --recal_file "$snp_out_1" --tranches_file "$snp_out_2" -o "$vqsr_out_1" ; then
echo "first step of variant recalibration: applying snp recalibration, was successfull"
else
echo "first step of variant recalibration: applying snp recalibration, failed!"
exit 1
fi

####STEP 2A: INDELS####
###################


indel_out_1=""$resultDir"/multisample.indel.model"
indel_out_2=""$resultDir"/multisample.indel.model.tranches"
indel_out_3=""$resultDir"/multisample.indel.model.plots.R"

declare -a vars=("$resultDir" "$java" "$gatk" "$multisample_in" "$known_indels" "$reference" "$known_snps" "$indel_out_1" "$indel_out_2" "$cpus" "$indel_out_3" "$vqsr_out_2" "$vqsr_out_1")
check_vars $vars "variant recalibraiton II" #function will abort script if one of those variables are unset
#known_indels = mills, known_snps = dbsnp 

if [ ! -f "$known_indels" ]; then echo "mills file could not be found, check path in config file" exit 1
fi

if
"$java" -jar "$gatk" -T VariantRecalibrator \
-R "$reference" --input "$multisample_in" --maxGaussians 4 -resource:mills,known=false,training=true,truth=true,prior=12.0 "$known_indels" -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "$known_snps" -an QD -an DP -an FS -an ReadPosRankSum -an MQRankSum --mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --recal_file "$indel_out_1" --tranches_file "$indel_out_2" -rscriptFile "$indel_out_3" -nt "$cpus"; then
echo "second step of variant recalibration: calculation of indel recalibration, was successfull"
else
echo "second step of variant recalibration: calculation of indel recalibration, failed!"
exit 1
fi

####STEP 2B: INDELS####
###################
#TODO maybe use ts_filter_level 90.0 to exclude more FP

if
"$java" -jar "$gatk" -T ApplyRecalibration -R "$reference" --input "$vqsr_out_1" --mode INDEL --ts_filter_level 99.0 --recal_file "$indel_out_1" --tranches_file "$indel_out_2" -o "$vqsr_out_2" ; then
echo "second step of variant recalibration: applying indel recalibration, was successfull"
else
echo "second step of variant recalibration: applying indel recalibration, failed!"
exit 1
fi

}

####################################################################################################
########################################HARD FILTERING DP AND AD ###############################
###################################################################################################
vcf_filter() {

echo "Filtering variants (first marking then selecting, excluding variants that didnt pass specified filter"

#filtering for depth (DP) and quality by depth (QD)
input_vcf="$1"
sample_ids="$2"
marked_vcf=""$resultDir"/multisample_marked_filter.vcf"

sample_ids_arr=(${sample_ids//:/ })
sample1="${sample_ids_arr[0]}"
sample2="${sample_ids_arr[1]}"
select_out=""$resultDir"/multisample_ready_variants_"$sample1"_"$sample2".vcf"

echo "Those are the samples you want to select after filtering:"
echo ""$sample1""
echo ""$sample2""

declare -a vars=("$resultDir" "$java" "$gatk" "$input_vcf" "$known_indels" "$reference" "$marked_vcf" "$select_out" "$sample1" "$sample2" ) 
check_vars $vars "variant selecting" #function will abort script if one of those variables are unset

echo "marking"
if "$java" -jar "$gatk" -T VariantFiltration -R "$reference" --variant "$input_vcf" --filterExpression "QD < 3.0" --filterName QDFilter --filterExpression "DP < 5" --filterName DPFilter -o "$marked_vcf"; then
echo " successfully marked input File according to the filter sepcified ( all variants with QD < 3 and DP < 5 ) are marked for filtering "
else
echo " VariantFiltration (marking) failed"
exit 1
fi

select_command_1=""$java" -jar "$gatk" -T SelectVariants -R "$reference" --variant "$marked_vcf" -o "$select_out" -selectType SNP --excludeNonVariants --excludeFiltered -sn "$sample1" -sn "$sample2" -select ""'""vc.getGenotype(\""
select_command_1b="$sample1"
select_command_1c="\").getAD().1 >= 2 && vc.getGenotype(\""
select_command_2="$sample1"
select_command_3="\").getDP() >= 5 && vc.getGenotype(\""
select_command_4="$sample2"
select_command_5="\").getAD().1 >= 2 && vc.getGenotype(\""
select_command_5b="$sample2"
select_command_5c="\").getDP() >= 5'"
select_command="$select_command_1""$select_command_1b""$select_command_1c""$select_command_2""$select_command_3""$select_command_4""$select_command_5""$select_command_5b""$select_command_5c"

echo "selecting"
if eval "$select_command" ; then
echo "selecting variants was successfull. Results can be found in "$select_out""
else
echo "selecting variants failed."
exit 1
fi



}
###########################################################
### COMMANDLINE PARSING ###
###########################################################



#looping through commandline args that are given by user

POSITIONAL=()
while [[ $# -gt 0 ]]
do
	key="$1"
	case $key in
		-h|--help)
		usage
		exit
		;;
		-o|--output)		
		resultDir="$2"
		if [ ! -d "$resultDir" ]; then
			echo "Specified output directory doen*t exist \n "
			exit 1
		else echo "output dir: "$resultDir""
		fi
		shift #shift arg
		shift #shift value
		;;
		-i|--input)
		Input_1="$2"
		echo "Your first Input file is: "$Input_1" "
		# somehow solution for case of only providing one inputfile ( e.g. for base_recal ()) 
		# in this case the next argument would starting with "-" 
		if [[ $3 == -* ]] || [[ $3 == --* ]] ; then 
			shift # shift argument
			shift # shift value 	
		else # if it doesnt start with "-" or "--", let*s assume it is another inputfile.
			Input_2="$3"
			echo " Your second Input file is: "$Input_2" "
			shift #shift arg
			shift #shift first value
			shift #shift second value
		fi
		;;
		-sID|--sampleID)
		sID="$2" #this is optional and if not provided by user, will be set with "set_vars()" later on
		shift
		shift
		;;
		 -VQSR_tr|--VQSR_threshold)
                VQSR_tr="$2" #this is optional and if not provided by user, will be set with "set_vars()" later on
                shift
                shift
		;;
		-c|--cpus)
		cpus="$2" #this is optional and if not provided by user, will be set with "set_vars()" later on
		shift
		shift
		;;
		-t|--tool)
		tool="$2"
		set_vars # check if optional parameters have been provided, otherwise setting variables to default values
		set_dependencies #parse config_wgs.xml and set all variables for dependencies 
		load_modules #loading tools with "module load <tool>"	
	#TODO function checking if all tools are working and no problems occur with file paths or verisons
		if [ $tool == "tr" ]; then 
			if parsing $Input_1 $Input_2; then #use parsing() to make sure vars are set and files exist
				if [ -z $resultDir ]; then 
					echo "you did not specify an output dir, a Result dir will be created where script was executed, trimmed reads will be saved in \"./Results\" "
					mkdir -p "$(pwd)/Results"
					resultDir="$(pwd)/Results"
				fi
			trimming $Input_1 $Input_2
			else
				echo "No files for trimming prodvided with -i (see -help)"
				exit 1
			fi
			exit
		elif [ $tool == "aln" ]; then
			if [ -z $Input_1 ] || [ -z $Input_2 ]; then
				printf "You did not specify Input files for alginment. \n Checking if trimmed reads are available in Result directory from previous trimming run:"
				Input_1=($resultDir/*_1P.fq*)
				Input_2=($resultDir/*_2P.fq*)
				if parsing "$Input_1" "$Input_2"; then
					echo "found trimmed reads in outpur dir \"Results\" "
				else
					echo "You did not specify any reads to align and no trimmed reads were found in Result dir, from previous trimming run"
					exit 1
				fi
			else
				if  ! parsing $Input_1 $Input_2; then
                                        echo "Provided Input files for alignment can*t be found or don*t exist!"
                                        exit 1
				fi
			fi
			alignment $Input_1 $Input_2
			exit
		elif [ "$2" == "brec" ]; then
			if [ -z $Input_1 ] || [ ! -f $Input_1 ]; then
				echo "No input file provided or cant be found! \n this is the path you provided: "$Input_1""
				exit 1
			else
				base_recal $Input_1
				exit
			fi
		elif [ "$2" == "varcall" ]; then		
			if [ -z $Input_1 ] || [ ! -f $Input_1 ]; then
				echo "No input file provided or cant be found! \n this is the path you provided: "$Input_1""
				exit 1	
			else
				varcall_discovery $Input_1
				exit
			fi
		elif [ "$2" == "varcall_compare" ]; then		
			if parsing "$Input_1" "$Input_2"; then
				varcall_given_alleles "$Input_1" "$Input_2"
				exit
			else
				echo "parsing input files failed"
				exit 1
			fi
		elif [ "$2" == "joint_call" ]; then	
			if [ -z $Input_1 ]; then echo "input string with vcf's to merge or g.vcf file to perform joint_call on, is not provided"
				exit 1
			fi
			# split input string by ":"
		        varArr=($(awk -F: '{$1=$1} 1' <<<"${Input_1}"))
			if [ ${#varArr[@]} -eq 1 ] ; then
				echo "You only provided one g.vcf file -> skipping merging"
				flag="0"
				joint_calling "$Input_1" "$flag"
				exit
			else	
				flag="1"
				echo "those are the files you want to merge and use for joint genotyping:"
				for i in "${varArr[@]}"
				do
				echo "$i"
				if [ ! -f "$i" ]; then
		        		echo "File "$i" doesn't exist, or cant be found"
		        		exit 1
				fi
				done
				joint_calling "$Input_1" "$flag"
			fi
			exit
		elif [ "$2" == "var_recal" ] ; then
			if [ -z "$Input_1" ] ; then 
				echo "no input provided"
				exit 1 
			fi
			if [ ! -f "$Input_1" ] ; then 
				echo "provided input file cannot be found" 
				exit 1
			 fi
			var_recal "$Input_1" "$VQSR_tr"
			exit 
		elif [ "$2" == "filter_var" ] ; then
			if [ -z "$Input_1" ] ; then 
				echo "no input provided"
				exit 1 
			fi
			if [ -z "$Input_2" ] ; then 
				echo "no sample ids specified to select variants from ( in the form of \"-i <multisample.vcf> <sampleid1_sampleid2>\")"
				exit 1 
			fi
			if [ ! -f "$Input_1" ] ; then 
				echo "provided input file cannot be found" 
				exit 1
			 fi
			vcf_filter "$Input_1" "$Input_2"
			exit
		else echo "Tool: "$tool" unknown. check your command \n " exit 1
		fi
		shift   
		shift
		;;
		--testing)
		select_command_1="java -jar gatk -T SelectVariants -R reference --variant markedvcf -o select_out --excludeNonVariants --excludeFiltered -sn sample1 -sn sample2 -select ""'""vc.getGenotype("
		select_command_1b="sample1"
		select_command_1c=").getAD().1 >= 2 && vc.getGenotype("
		select_command_2="sample1"
		select_command_3=").getDP() >= 5 && vc.getGenotype("
		select_command_4="sample2"
		select_command_5=").getAD().1 >= 2 && vc.getGenotype("
		select_command_5b="sample2"
		select_command_5c=").getDP() >= 5 '"
	select_command="$select_command_1""$select_command_1b""$select_command_1c""$select_command_2""$select_command_3""$select_command_4""$select_command_5""$select_command_5b""$select_command_5c"
		echo "$select_command"
                #set_dependencies
		#load_modules
		exit 
                shift
                ;;
		*) #unknown option
		echo "unknown option: "$1"\nRun ${0##*/} -h for help.">&2
		exit 1
	#	POSITIONAL+=("$1") # save in array
	#	shift #past argument
	#	;;
	esac
done

set -- "${POSITIONAL[@]}" #restore positional parameters	

