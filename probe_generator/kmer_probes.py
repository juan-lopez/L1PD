#!/usr/bin/env python3
"""This script processes SAM files to find the k-mers that have alignments
in the indicated ORF but the least amount of alignments outside of it.
These k-mers will be candidates to use as probes.
"""

import argparse
import json
import os
import re
import sys

# Import shared L1Base2 module
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'L1PD_files'))
from merge.algorithm import kmer_frequency, heap_merge_kmers
import l1base2

from Bio import SeqIO

threshold = None # Script argument; Should use same threshold when filtering
# We will assume a file name format so that we may extract the corresponding chromosome
# from the file name, and then include it in the results.
chromosome = None
DEBUG = False
count = 1
KMER_FREQ_STEP = 50


def load_file_data(SAMFileName):
	"""Return nested dictionary with alignment pos in SAM file, and k-mer size.

	First key of dictionary is the chromosome, second key is the k-mer name.
	"""
	kmerSize = -1 # K-mer size will be extracted from CIGAR in SAM file lines
	with open(SAMFileName) as fh:
		kmerDict = dict() # Dictionary of dictionaries ([chromosome][kmer])
		for line in fh:
			if line.startswith("@"): # Skip header lines
				continue
			pieces = line.split("\t") # SAM files have tab-separated columns
			kmerName = pieces[0]
			chrm = pieces[2] # chromosome
			cigar = pieces[5] # CIGAR string
			if kmerSize == -1: # If haven't detected k-mer size yet
				# Don't want to do much math, simply find a line with all M's
				# Note that we're assuming all k-mers have the same size.
				m = re.match(r'(\d+)M$', cigar)
				if m:
					kmerSize = int(m.group(1))
			if chrm not in kmerDict:
				kmerDict[chrm] = dict()
			# If the FLAG field has 5th bit (16) on, seq is reverse complemented.
			# In that case we store the position as negative to indicate so.
			kmerPos = int(pieces[3]) if not int(pieces[1])&16 else -int(pieces[3])
			# It's possible for a position to appear more than once for a chromosome,
			# so only add it once to the list of positions.
			if kmerName not in kmerDict[chrm]: # First pos for this k-mer in this chrm
				kmerDict[chrm][kmerName] = [kmerPos]
			elif kmerPos not in kmerDict[chrm][kmerName]:
				kmerDict[chrm][kmerName].append(kmerPos)
	return (kmerDict, kmerSize)


def discard_not_in_ORF(kDict, Prefix, Verbose):
	"""Remove from dict all k-mers that don't have any hits inside the appropriate ORF."""
	L1sByChar = l1base2.get_L1s_by_chr()
	# We want to keep track of which L1s' ORFs were matched by which k-mers
	ORFsMatched = dict() # Another nested dict with L1 numbers
	ORFsMatched["All"] = dict() # Store stats for all chromosomes
	for chrm in kDict: # For every chromosome found in SAM file
		if chrm not in L1sByChar: # No L1s for this chromosome
			continue
		kmersToRemove = list()
		for kmer in kDict[chrm]: # For every k-mer alignment in that chromosome
			# Compare if any pos for k-mer is within ORF of any L1 for this chrm
			posKmerFound = list() # Pos in this L1's ORF where k-mer was found
			for L1Num in L1sByChar[chrm]: # Compare with every L1 for that chrm
				foundInThisL1 = False
				for pos in kDict[chrm][kmer]: # For every pos in k-mer's alignments
					if pos_inside_ORF(pos, L1Num, Prefix):
						posKmerFound.append( (Prefix, chrm, kmer, L1Num, pos) )
						if foundInThisL1:
							if Verbose:
								print("FOUND MORE THAN ONCE IN ORF OF SAME L1:",posKmerFound,file=sys.stderr)
							break # No need to process this k-mer any more
						else:
							# Keep track of which L1's ORF was matched, since we will later filter
							# out k-mers based on which chrms and/or L1's were not matched.
							if chrm not in ORFsMatched:
								ORFsMatched[chrm] = dict()
							ORFsMatched[chrm][kmer] = ORFsMatched[chrm].get(kmer,list())+[L1Num+1]
							ORFsMatched["All"][kmer] = ORFsMatched["All"].get(kmer,list())+[L1Num+1]
						foundInThisL1 = True
			if len(posKmerFound) == 0: # No positions for this k-mer matched, so discard it
				kmersToRemove.append(kmer) # Can't remove it mid-iteration
		for kmer in kmersToRemove:
			del kDict[chrm][kmer]
	return ORFsMatched


#def discard_outside_ranges(kDict, k, ranges):
#	"""Remove from dict all k-mer hits outside ranges with >= 95% col identity."""
#	L1sByChar = l1base2.get_L1s_by_chr()
#	for chrm in kDict: # For every chromosome found in SAM file
#		if chrm not in L1sByChar: # No L1s for this chromosome
#			continue
#		kmersToRemove = list()
#		for kmer in kDict[chrm]: # For every k-mer alignment in that chromosome
#			foundPosInRange = False
#			for pos in kDict[chrm][kmer]:
#				if pos_inside_ranges(pos, k, ranges[chrm]):
#					foundPosInRange = True
#					break
#			if not foundPosInRange: # K-mer has no hits within the ranges
#				kmersToRemove.append(kmer) # Can't remove it mid-iteration
#		for kmer in kmersToRemove: # Remove k-mers with no hits in ranges
#			del kDict[chrm][kmer]


def discard_lower_L1_coverage(kDict, ORFsMatched):
	"""Remove k-mers that don't have hits across ALL L1s."""
	kmersToRemove = { kmer for kmer in ORFsMatched["All"] if len(ORFsMatched["All"][kmer]) < len(l1base2.get_CSV()) }

	for chrm in kDict:
		# Remove the k-mers from kmersToRemove that are in this chromosome
		for kmer in kmersToRemove.intersection(set(kDict[chrm].keys())):
			del kDict[chrm][kmer]


def pos_inside_ORF(pos, L1Num, Prefix):
	"""Determine whether a position lies within the ORF of a particular L1."""
	CSV = l1base2.get_CSV()
	ALU = Prefix # Look within a specific ORF
	if (pos >= 0 and CSV[L1Num]["Strand"] == "1" and # Forward strand
			 int(CSV[L1Num]["Start"])+int(CSV[L1Num][ALU+" Start"]) <= pos and
			 int(CSV[L1Num]["Start"])+int(CSV[L1Num][ALU+" End"])   >= pos) or \
			(pos < 0 and CSV[L1Num]["Strand"] == "-1" and # Reverse strand
			 int(CSV[L1Num]["Start"])-int(CSV[L1Num][ALU+" Start"]) >= -pos and
			 int(CSV[L1Num]["Start"])-int(CSV[L1Num][ALU+" End"])   <= -pos):
		return True
	else:
		return False


#def pos_inside_ranges(pos, k, ranges):
#	"""Determines whether the entire k-mer lies within any of the ranges."""
#	# Ranges are sorted, so no need to check them all
#	rangeIndex = 0
#	if pos >= 0: # Forward strand
#		while rangeIndex < len(ranges) and ranges[rangeIndex][0] <= pos:
#			if pos+k-1 <= ranges[rangeIndex][1]: # Loop already tests lower end
#				return True
#			rangeIndex += 1
#	else: # Reverse strand
#		pos = -pos
#		while rangeIndex < len(ranges) and ranges[rangeIndex][0] <= pos-k+1:
#			if pos <= ranges[rangeIndex][1]: # Loop already tests lower end
#				return True
#			rangeIndex += 1
#	return False


def print_min_spread(kPosDict,ORFsMatched, KmerFile, AlignedKmersFile, k, Prefix, Verbose, CSVFile, Merge, is_line_1, max_ambiguity_treshold):
	kPosDict["All"] = dict() # Used to store total amount of pos for all chrm
	tupList = list()
	append = tupList.append # Local variable; hopefully more efficient

	for chrm in kPosDict:
		if chrm == "All":
			continue

		if Verbose:
			print("Chromosome", chrm,file=sys.stderr)

		if chrm in ORFsMatched or not is_line_1:
			tupList.clear()
			for kmer in kPosDict[chrm]:
				if is_line_1:
					append( (len(kPosDict[chrm][kmer]), kmer, ORFsMatched[chrm][kmer]) ) 
				
				else:
					append( (len(kPosDict[chrm][kmer]), kmer) )
				
				kPosDict["All"][kmer] = kPosDict["All"].get(kmer,0) + len(kPosDict[chrm][kmer])

			tupList.sort()
			if len(tupList) == 0: # This should not happen!
				print("tupList is empty! kPosDict["+chrm+"] =", kPosDict[chrm])

			if Verbose:
				print(*tupList, sep='\n', file=sys.stderr)

		else: # Some chromosomes don't have L1s
			if Verbose:
				print("No ORFs matched.  No L1s in this chromosome?",file=sys.stderr)

	# Now print in descending order of ORFs matched, and ascending order
	# of "spread", across all chromosomes.
	# At the same time, find if k-mer is present in aligned version.
	if Verbose:
		print("Global",file=sys.stderr)

	aligned_kmers = list(SeqIO.parse(AlignedKmersFile,"fasta"))
	# The SAM file had the k-mer names, but now we need the actual k-mers
	rd = SeqIO.to_dict(SeqIO.parse(KmerFile,"fasta"))
	tupList = list()

	if is_line_1:
		for kmer in kPosDict["All"]:
			# Try to find this k-mer in one of the ORFs
			kmerPosInORF = -1
			ORFIndex = 0
			while ORFIndex < len(aligned_kmers) and kmerPosInORF == -1:
				kmerPosInORF = aligned_kmers[ORFIndex].seq.find(rd[kmer].seq)
				ORFIndex += 1
			tupList.append( (len(ORFsMatched["All"][kmer]), kPosDict["All"][kmer], kmer, kmerPosInORF) )
		tupList.sort(reverse=True, key=lambda x: (x[0], -x[1]))

	else:
		# Iterates over kmer names
		for kmer in kPosDict["All"]:
			if kPosDict["All"][kmer] > 0:
				# Try to find this k-mer in one of the aligned alu sequences
				kmerPos = -1
				AluIndex = 0
				while AluIndex < len(aligned_kmers) and kmerPos == -1:
					kmerPos = aligned_kmers[AluIndex].seq.find(rd[kmer].seq)
					AluIndex += 1

				if kmerPos != -1:
					tupList.append( (kPosDict["All"][kmer], kmer, kmerPos) )
		tupList.sort(reverse=True, key=lambda x: (x[0], -x[2]))
	finalTupList = list()

	if Merge:
		finalTupList = heap_merge(tupList, k, max_ambiguity_treshold, Prefix)
	
	# Ensure only non-overlapping k-mers are used
	else:
		for tup in tupList:
			# Check last component of tuple to see if k-mer was found in Fixed file
			tupPos = tup[-1]
			if tupPos > -1: # Component has result of find method
				# Check if this position doesn't conflict with previous positions
				conflict = False
				for finalTup in finalTupList:
					finalTupPos = finalTup[-1]
					if finalTupPos <= tupPos + k - 1 and tupPos <= finalTupPos + k - 1: 
						conflict = True
						break
				if not conflict:
					if is_line_1:
						tup = tup[1:] # Omit the L1 count
					finalTupList.append(tup) 
	if Verbose:
		print(*finalTupList, sep='\n', file=sys.stderr)

	# Finally, print these final k-mers in FASTA format
	count = 1
	# File containing default consensus (identity percentage 0) for merge
	all_bases = ""
	with open("all_bases.txt", "r") as file:
		all_bases = file.readline().strip()

	print(all_bases)

	for tup in sorted(finalTupList,key=lambda x: x[2]):
		#if tup[1] not in rd: TODO DEBUG
		#	print("tup:", tup)
		#	print("The following key is not in rd:", tup[1])

		# extracting the actual bases
		if Merge:
			# from the default consensus (identity percentage 0) for merge
			sequence = all_bases[tup[2]:tup[2] + k]

		else:
			# from Kmer name when not merging
			sequence = rd[tup[1]].seq

		if Prefix == "ORF1":
			print(f'>{Prefix}_{k}mers_{count}',str(tup[2])+'\n'+sequence)
		elif Prefix == "ORF2": # if ORF2, add the Avg ORF1 size and Inter-ORF size
			orf1_inter_sum_result = sum_orf1_inter(CSVFile)
			print(f'>{Prefix}_{k}mers_{count}',str(tup[2]+orf1_inter_sum_result)+'\n'+sequence)
		else:
			#Output format: >kmer-name kmer position \n kmer
			print(f'>{Prefix}_{k}mers_{count}',str(tup[2])+'\n'+sequence)
		count += 1


def sum_orf1_inter(CSVFile):
	"""Adds the lenght of the ORF1 and Inter elements
	
	Parameters:
	CSVFile - Contains the path where the .csv files that have the meta data are

	Return:
		integer number resulting from the sum of ORF1 and Inter lengt.
	
	"""
	l1base2.load_CSV_file(CSVFile)
	avgLen = l1base2.AvgLengths
	sum = int(avgLen["ORF1"]) + int(avgLen["Inter"])
	return sum


def get_count_kmer_algnmnt(kDict):
	kmerCount = alignmentCount = 0
	kmerSet = set()
	for chrm in kDict:
		for kmer in kDict[chrm]:
			kmerSet.add(kmer)
			alignmentCount += len(kDict[chrm][kmer])
	return (len(kmerSet), alignmentCount)


def merge_overlap(tupList, k):
    """
    Merge ONLY kmers that overlap (no gaps allowed).
    
    Parameters:
        tupList: List of tuples (L1Count, hits, kmer_name, position)
        k: Original k-mer size (e.g., 50)
    
    Returns:
        List of tuples (start_position, merged_length) for merged regions
    """
    if not tupList:
        return []
    
    # Sort by position (4th element of tuple)
    tupList.sort(key=lambda x: x[3])
    
    merged = []
    current_start = tupList[0][3]
    current_end = current_start + k  # Each kmer is length k
    
    for i in range(1, len(tupList)):
        next_pos = tupList[i][3]
        next_end = next_pos + k
        
        # Check if there's ANY overlap (next kmer starts before current ends)
        # Overlap exists if: next_pos < current_end
        # (Use <= if you want to merge adjacent kmers with 0 gap)
        if next_pos <= current_end:
            # They overlap! Extend current region
            current_end = max(current_end, next_end)
        else:
            # No overlap - save current merged region
            merged_length = current_end - current_start
            if merged_length >= k:  # Keep only if at least original k size
                merged.append((current_start, merged_length))
            
            # Start new region
            current_start = next_pos
            current_end = next_end
    
    # Don't forget the last region!
    final_length = current_end - current_start
    if final_length >= k:
        merged.append((current_start, final_length))
    
    return merged


def find_gaps_between_merged_kmers(merged_tup_list):
    """
    Calculate the gaps (non-kmer regions) between merged kmers.
    
    Parameters:
        merged_tup_list: List of (start_position, length) tuples for merged kmers
    
    Returns:
        List of (gap_start, gap_length) tuples for gaps between kmers
    """
    if len(merged_tup_list) < 2:
        return []  # Need at least 2 kmers to have gaps between them
    
    non_kmers = []
    
    for i in range(len(merged_tup_list) - 1):
        start_a, len_a = merged_tup_list[i]
        start_b, len_b = merged_tup_list[i + 1]
        
        # End position of first kmer
        end_a = start_a + len_a
        
        # Gap starts right after end_a, ends right before start_b
        gap_start = end_a
        gap_end = start_b
        
        # Calculate gap length (only if there's actually a gap)
        gap_length = gap_end - gap_start
        
        if gap_length > 0:
            non_kmers.append((gap_start, gap_length))
        # If gap_length <= 0, kmers overlap or touch - no gap
    
    return non_kmers

def split_kmers(tupList, newKmerSize, orf):
    """
    Cuts the tup list in to individual k-mers of length: newKmerSize
    """
    splitKmers = []
    for start, length in tupList:
        while length - newKmerSize > 0:
            kmerName = f"{orf}_{newKmerSize}mers_"
            splitKmers.append((newKmerSize, kmerName, start))
            length -= newKmerSize
            start += newKmerSize

    return splitKmers

def heap_merge(tupList, k, max_ambiguity_treshold, prefix):
	max_ambiguity_treshold /= 10
	merged_tup_list = merge_overlap(tupList, k)
	non_kmers = find_gaps_between_merged_kmers(merged_tup_list)
	final_kmers = heap_merge_kmers(merged_tup_list, non_kmers, max_ambiguity_treshold)
     
	print("------------------------------------------------")
	print("sliding window merged kmers", len(merged_tup_list))
	# Filter by k-mer threshold
	final_kmers = [(start, length) for start, length in final_kmers if length >= k]
	print("Final kmers")
	print(final_kmers)
	'''
	print("------------------------------------------------")
	print("Heap merge kmers", len(final_kmers))
	print("Min kmer length")
	print(min(final_kmers, key=lambda x: x[1])[1])
	print("Max kmer length")
	print(max(final_kmers, key=lambda x: x[1])[1])
	print("Average length")
	print(sum([kmer[1] for kmer in final_kmers]) / len(final_kmers))
	print(kmer_frequency(final_kmers, 50))
	final_kmers = split_kmers(final_kmers, k, prefix)
	print("Heap merge kmers after split") 
	'''
	final_kmers = split_kmers(final_kmers, k, prefix)
	print("Final kmers after split")
	print(len(final_kmers))
	print(final_kmers)

	return final_kmers


# Function below is obsolete, but left for future reference for JSON storage
# This way the data is loaded from JSON files and not generated again.
#def main_new(SAMFile, KmerFile, AlignedORFsFile, CSVFile):
#	"""OBSOLETE Run the algorithm, but load the data from JSON files."""
#	ORFNum = int(KmerFile[3]) # KmerFile must start with ORF1 or ORF2
#	with open('kmerPosDict_ORF'+str(ORFNum)+'.json') as fp:
#		kmerPosDict = json.load(fp)
#	with open('ORF'+str(ORFNum)+'sMatched.json') as fp:
#		ORFsMatched = json.load(fp)
#	kmerSize = 50
#	print_min_spread(kmerPosDict, ORFsMatched, KmerFile, AlignedORFsFile, kmerSize)


def main(SAMFile, KmerFile, AlignedORFsFile, CSVFile, Prefix, Verbose, Merge, is_line_1, max_ambiguity_treshold):
	if is_line_1:
		l1base2.load_CSV_file(CSVFile)
	kmerDict = SeqIO.to_dict(SeqIO.parse(KmerFile,"fasta"))
	kmer_count = len(kmerDict)

	# Generate a dictionary with alignments from SAM file.
	kmerPosDict, kmerSize = load_file_data(SAMFile)
	if Verbose:
		if kmerSize == -1: # Unable to detect k-mer size
			print("Unable to detect k-mer size in", SAMFile, "; skipping file.",file=sys.stderr)
			exit()

	oldKCount, oldACount = get_count_kmer_algnmnt(kmerPosDict)
	if Verbose:
		print(str(oldKCount)+"/"+str(kmer_count),"("+str(100*oldKCount/kmer_count)+"%) k-mers with", oldACount, "alignments",file=sys.stderr)
	
	ORFsMatched = 0
	if is_line_1:
		if oldKCount > 0:
			# Discard k-mers without alignments inside targeted ORF
			ORFsMatched = discard_not_in_ORF(kmerPosDict, Prefix, Verbose)
			newKCount, newACount = get_count_kmer_algnmnt(kmerPosDict)
			if Verbose:
				print("After discarding ^In, {}/{} ({}%) k-mers left with {}/{} ({}%) alignments".format(newKCount,oldKCount,100*newKCount/oldKCount,newACount,oldACount,100*newACount/oldACount),file=sys.stderr)

			if newKCount > 0:
				discard_lower_L1_coverage(kmerPosDict, ORFsMatched)
				oldKCount, oldACount = newKCount, newACount
				newKCount, newACount = get_count_kmer_algnmnt(kmerPosDict)
				if Verbose:
					print("After discarding lower coverage, {}/{} ({}%) k-mers left with {}/{} ({}%) alignments".format(newKCount,oldKCount,100*newKCount/oldKCount,newACount,oldACount,100*newACount/oldACount),file=sys.stderr)
	else:
		ORFsMatched = kmerDict
	## Print k-mers in ascending order of amount of alignments
	# We can store the k-mer alignments in JSON files to avoid having
	# to repeat these steps.  In that case, we could then load the
	# JSON files and proceed with the rest of the process.
	#with open('kmerPosDict_'+Prefix+'.json','w') as fp:
		#json.dump(kmerPosDict,fp,indent=3)
	#with open(Prefix+'sMatched.json','w') as fp:
		#json.dump(ORFsMatched,fp,indent=3)
	print_min_spread(kmerPosDict, ORFsMatched, KmerFile, AlignedORFsFile, kmerSize, Prefix, Verbose, CSVFile, Merge, is_line_1, max_ambiguity_treshold)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Find k-mers with aligment in ORF but least alignments elsewhere")
	parser.add_argument("SAM", help="Full path to SAM file with alignments")
	parser.add_argument("KmerFile", help="Full path to FASTA file with all k-mers in ORF")
	parser.add_argument("AlignedORFsFile", help="Full path to FASTA file with aligned ORFS (all having same length)")
	parser.add_argument("L1BaseCSV", help="Full path to L1Base CSV file with L1 data")
	parser.add_argument("Prefix", help="Prefix for the files generated")
	parser.add_argument('-v', '--verbose', help="Add verbosity so that output contains additional information sent to the standard error output", action='store_true')
	parser.add_argument('-m', help="Enables the merging of overlapping kmers which are sent to the standard error output",action='store_true')
	parser.add_argument("-q", type=int, default=1, help="Type of sequence: LINE or Other")
	parser.add_argument("-p", type=int, default=1, help="Max ambiguity treshold")
	args = parser.parse_args()
	main(args.SAM, args.KmerFile, args.AlignedORFsFile, args.L1BaseCSV, args.Prefix, args.verbose, args.m, args.q, args.p)