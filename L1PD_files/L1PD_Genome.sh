#!/bin/sh -

########################################################################
# L1PD_Genome.sh
#
# Execute the L1PD algorithm on a genome, specified as argument.
#
# Author: Juan O. Lopez (juano.lopez@upr.edu)
# License: Creative Commons Attribution-ShareAlike 4.0
# http://creativecommons.org/licenses/by-sa/4.0/legalcode
########################################################################

# "Constants"
PROBESFASTA="L1PD_50mers.fasta"
#REFERENCECSV="/anakena/1000Genome/data/GRCh38DH.csv"
#SCALEFACTOR=100

# Initialize variables to default values
edit_distance=20
dist_threshold=700
min_amt_kmers=9

while getopts p:e:d:m:f:h opt ; do
	case $opt in
		p) prefix=$OPTARG ;;
		e) edit_distance=$OPTARG ;;
		d) dist_threshold=$OPTARG ;;
		m) min_amt_kmers=$OPTARG ;;
		f) file=$OPTARG ;;
		h) usage
			exit 0
			;;
	esac
done

if [ -z $file ] ; then
	echo "ERROR: Genome in FASTA format is missing!"
	exit 1
fi
if [ -z $prefix ] ; then
	# If no prefix was specified, use the filename prefix (up to first
	# period) as the filename (with different extensions) for the
	# resulting files.  If prefix has a space, only first word is used.
	# Note: Use quotes around the file name to allow spaces.
	algnBase=$(basename $file)
	prefix=${algnBase%%.*}
fi

# Testing
echo Testing... L1PD_Genome -f $file -e $edit_distance -t $dist_threshold -p $prefix -m $min_amt_kmers

dir_path=$(dirname $0)
# If .index file doesn't exist, then index the genome
if [ ! -f $file.index ] ; then
	echo "Index file doesn't exist; indexing genome..."
	mrfast --index $file
fi
mrfast --search $file --seq $dir_path/$PROBESFASTA -e $edit_distance -o probes_$prefix.sam
python3 $dir_path/L1PD.py probes_$prefix.sam $dir_path/$PROBESFASTA -t $dist_threshold -m $min_amt_kmers > $prefix.gff3
#rm probes_$prefix.sam

# Calculate CNV
#PCSUB=$(wc -l ${prefix}.gff3 | sed 's/\s.*//')
#PCREF=$(wc -l $REFERENCECSV | sed 's/\s.*//')
#CNV=$(echo "($PCSUB - $PCREF)*$SCALEFACTOR/$PCREF" | bc -l) 
#echo "Subject pattern count: $PCSUB"
#echo "Reference pattern count: $PCREF"
#echo "CNV: $CNV"
