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

# We use a function to download the FASTQ if necessary.
# Additionally, it what was downloaded was a .gz file, we also unzip it.
retrieve_FASTQ()
{
	# If the filename starts with ftp:// then we download it.
	FASTQ=$1
	if beginswith "ftp://" "$FASTQ" ; then
		wget $FASTQ
		# Now that it's downloaded, strip off URL and keep only filename
		FASTQ=$(basename $FASTQ)
		if endswith ".gz" "$FASTQ" ; then
			# First gunzip and then remove '.gz' from the filename
			gunzip $FASTQ
			FASTQ=${FASTQ:0:-3}
		fi
	fi
}

## Retrieve both FASTQs and update the names to the local filenames
#retrieve_FASTQ $file1
#file1=$FASTQ
#retrieve_FASTQ $file2
#file2=$FASTQ

# If .bwt file doesn't exist, then index the reference genome
if [ ! -f $reference.bwt ] ; then
	echo "BWT file doesn't exist; indexing reference genome..."
	bwa index -a bwtsw $reference
fi

SAM=$prefix.sam
bwa mem -t $threads $reference $file1 $file2 > $SAM

# Create BAM from SAM, sort, fix mate information and add the MD tag
BAM=$prefix.bam
cmd="samtools fixmate -u -O bam -@$threads $SAM -"
cmd="$cmd | samtools sort -u -@$threads -T $prefix -"
cmd="$cmd | samtools calmd -b -@$threads - $reference > $BAM"
$cmd
# SAM file takes up too much space
rm $SAM

# Now L1PD_BAM_CRAM.sh takes over.
