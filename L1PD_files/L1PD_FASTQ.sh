#!/bin/sh -

########################################################################
# L1PD_FASTQ.sh
#
# This script aligns the paired-end FASTQ files to the reference genome
# using BWA and then converts the SAM file to a BAM file.
#
# Author: Juan O. Lopez (juano.lopez@upr.edu)
# License: Creative Commons Attribution-ShareAlike 4.0
# http://creativecommons.org/licenses/by-sa/4.0/legalcode
########################################################################

usage() {
	cat << EOF
Usage: L1PD_FASTQ.sh -f file1.fq file2.fq -r ref.fa [-t threads] [-p prefix] [-h]

	-f  Paired-end FASTQ files [REQUIRED]
	-r  Reference genome in FASTA format [REQUIRED]
	-t  Number of additional threads
	    [default: 0]
	-e  Max. edit distance
	    [default: 20]
	-d  Distance difference threshold for establishing patterns
	    [default: 700]
	-p  Prefix text to use for all file names generated
	    [default: Filename prefix of FASTQ/BAM/CRAM file]

	-h  Usage
EOF
}

# https://stackoverflow.com/questions/2172352/in-bash-how-can-i-check-if-a-string-begins-with-some-value#18558871
beginswith() { case $2 in "$1"*) true;; *) false;; esac; }
endswith() { case $2 in *"$1") true;; *) false;; esac; }

while getopts hr:p:t:f: opt ; do
	case $opt in
		r) reference=$OPTARG ;;
		p) prefix=$OPTARG ;;
		t) threads=$OPTARG ;;
		h) usage
			exit 0
			;;
		f) 
			## $OPTIND has the index of the _next_ parameter
			#eval "file1=\${$OPTIND}"
			## "shift" getopts' index
			#OPTIND=$((OPTIND+1))
			file1=$OPTARG
			# Script can accept a single alignment file (BAM/CRAM), or two
			# FASTQ files, so test whether there's an additional file argument.
			if [ $OPTIND -le $# ] && case $(eval "echo \${$OPTIND}") in -*) false;; *) true; esac; then
				eval "file2=\${$OPTIND}"
				OPTIND=$((OPTIND+1))
			else
				usage
				exit 1
			fi
			;;
	esac
done

if [ -z $file1 ] ; then
	# No filenames were provided
	usage
	exit 1
fi
if [ -z $prefix ] ; then
	# If no prefix was specified, use the filename prefix (up to first
	# period) as the filename (with different extensions) for the
	# resulting files.  If prefix has a space, only first word is used.
	# Note: Use quotes around the file name to allow spaces.
	algnBase=$(basename $file1)
	prefix=${algnBase%%.*}
fi

# If threads argument not received, use 0 (bcftools' default)
: ${threads:=0}

# Testing
echo Testing... L1PD_FASTQ -f $file1 $file2 -r $reference -t $threads -p $prefix

# We use a function to download the FASTQ if necessary.
# Additionally, it what was downloaded was a .gz file, we also unzip it.
retrieve_FASTQ()
{
	# If the filename starts with ftp:// then we download it.
	FASTQ=$1
	#if [[ $FASTQ = 'ftp://'* ]] ; then
	if beginswith "ftp://" "$FASTQ" ; then
		echo "Downloading $FASTQ with wget..."
		wget $FASTQ
		# Now that it's downloaded, strip off URL and keep only filename
		FASTQ=$(basename $FASTQ)
		#if [[ $FASTQ = *'.gz' ]] ; then
		if endswith ".gz" "$FASTQ" ; then
			# First gunzip and then remove '.gz' from the filename
			echo "Gunzipping $FASTQ..."
			time gunzip $FASTQ
			FASTQ=${FASTQ:0:-3}
		fi
	fi
}

# Retrieve both FASTQs and update the names to the local filenames
retrieve_FASTQ $file1
file1=$FASTQ
retrieve_FASTQ $file2
file2=$FASTQ

# If .bwt file doesn't exist, then index the reference genome
if [ ! -f $reference.bwt ] ; then
	echo "BWT file doesn't exist; indexing reference genome..."
	time bwa index -a bwtsw $reference
fi

# Not using Backtrack; using MEM instead
## Generate Suffix Array Indexes (SAI) files
## q is the parameter for read trimming; 1000 Genomes used 15, so do the same
#SAI1=$(basename $file1).sai
#SAI2=$(basename $file2).sai
#echo "BWA ALN #1"
#time bwa aln -q 15 -t $threads $reference $file1 > $SAI1
#echo "BWA ALN #2"
#time bwa aln -q 15 -t $threads $reference $file2 > $SAI2

# Generate SAM files
#echo "BWA SAMPE"
SAM=$prefix.sam
#time bwa sampe $reference $SAI1 $SAI2 $file1 $file2 > $SAM
echo "BWA MEM"
time bwa mem -t $threads $reference $file1 $file2 > $SAM

# Create BAM from SAM, sort, fix mate information and add the MD tag
#echo "SAM -> BAM, sort, fix mate, sort, add MD"
echo "fix mate, sort, add MD"
BAM=$prefix.bam
# Updated according to https://www.htslib.org/workflow/
#Old 1000Genome style: time ( samtools view -bSu $SAM | samtools sort -n -o - - | samtools fixmate - - | samtools sort -o - - | samtools calmd -b --threads $threads - $reference > $BAM )
#Only works in bash: set -o pipefail
cmd="samtools fixmate -u -O bam -@$threads $SAM -"
cmd="$cmd | samtools sort -u -@$threads -T $prefix -"
cmd="$cmd | samtools calmd -b -@$threads - $reference > $BAM"
time sh -c "$cmd"
# SAM file takes up too much space
rm $SAM
#time ( samtools fixmate -O bam -@$threads $SAM - | \
	    #samtools sort -u -@$threads -T $prefix - | \
		 #samtools calmd -b -@$threads - $reference > $BAM )

## Execute PC_BAM_CRAM.sh in the same directory as this script
#script_path=$(dirname "$0")
#$script_path/PC_BAM_CRAM.sh $BAM $reference $threads
