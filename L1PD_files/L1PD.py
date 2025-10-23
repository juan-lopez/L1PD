#!/usr/bin/env python3
"""This class aims to find patterns of hits of our k-mer probes that
have particular distances between them, providing possible L1s.
The algorithm is named LINE-1 Pattern Detection (L1PD).
The patterns that meet the established criteria are considered to be L1s,
and they are output into a file in GFF3 format.
NOTE: Only chromosomes of the form chrm#, where # is 1-22, are considered.
"""

import argparse
import bisect
import re
import time
import os
import sys
import dbRIP
import l1base2

threshold = None # Script argument; Should use same threshold when filtering
DEBUG = None # Script argument for printing debugging messages
CSVOutput = None # Boolean to determine whether to output in GFF3 or CSV format


def same_strand(pos1, pos2):
    """Determine if both positions are on the forward or reverse strand."""
    # Reverse complemented sequences have their position as negative numbers,
    # so use exclusive or (XOR) to compare if they both have the same sign.
    return (pos1 ^ pos2) >= 0


def find_L1s(kmerPosDict, kmerNameList, distList, minKmers, CSVdir, sequence_type):
    """Find patterns of hits that constitute LINE-1s.

    Parameters:
        kmerPosDict  - Contains the positions of the k-mer alignments
        kmerNameList - Names of the k-mers extracted from SAM; keys for kmerPosDict
        distList     - Positions of k-mers (to a fixed reference point)
        minKmers     - Minium amount of k-mers that must be found for a LINE-1

    Return:
      Dictionary with chromosomes as keys and L1 count as values
    """
    L1sPerChrmCount = dict() # Store the L1 count per chromosome

    amtKmers = len(kmerNameList)
    if not CSVOutput:
        print("##gff-version 3.1.25") # Copied from GFF3 Spec page

    for chrm in kmerPosDict:
        L1sPerChrmCount[chrm] = 0
        # Since we don't need to have all k-mers to have found a LINE-1, we need to be
        # careful that we don't use positions more than once when trying to generate
        # LINE-1s that don't include the first k-mer.
        usedKmerPos = [ set() for _ in range(len(kmerNameList)) ]
        if DEBUG:
            print("chrm =", chrm)
        # We want to get as big a LINE-1 as possible, so we start with the first k-mer
        # and keep trying different start k-mers, but making sure we don't repeat
        # positions (i.e. don't consider a segment of a previously recorded LINE-1).
        for startKmerIndex in range(0, amtKmers-minKmers+1):
            startKmerName = kmerNameList[startKmerIndex]
            if DEBUG:
                print("startKmerIndex =", startKmerIndex," startKmerName =", startKmerName)
            if startKmerName not in kmerPosDict[chrm]:
                if DEBUG:
                    print("startKmerName not in kmerPosDict["+chrm+"]")
                continue
            # Our start k-mer could have multiple positions depending on the mapping results
            # Try each one to see if we detect a LINE-1
            for startKmerPos in kmerPosDict[chrm][startKmerName]:
                if DEBUG: print("startKmerPos =", startKmerPos)
                if startKmerPos in usedKmerPos[startKmerIndex]: # Used in prev L1?
                    if DEBUG: print("Already used!")
                    continue
                # Try to find all possible LINE-1s beginning with the start k-mer
                L1 = [ (startKmerIndex, startKmerPos) ] # List of tuples
                for kmerIndex in range(startKmerIndex + 1, amtKmers):
                    if kmerNameList[kmerIndex] not in kmerPosDict[chrm]: # No hits for this k-mer
                        continue
                    # Calculate expected position of this k-mer so that it maintains the
                    # same distance to the start k-mer as in the FASTA file. 
                    distBetKmers = distList[kmerIndex] - distList[startKmerIndex] # Dist between k-mers
                    expectedPos = startKmerPos + distBetKmers # Works for forward and reverse
                    if DEBUG:
                        print("kmerIndex =", kmerIndex, " distBetKmers =", distBetKmers," expectedPos =", expectedPos)
                    # Verify whether this position has already been used for this k-mer in another LINE-1
                    if expectedPos not in usedKmerPos[kmerIndex]:
                        tempList = kmerPosDict[chrm][kmerNameList[kmerIndex]]
                        # bisect module does binary search faster!
                        actualPosIndex = bisect.bisect_left(tempList, expectedPos - threshold)
                        # Check whether the position is within the threshold of the expected position
                        if actualPosIndex != len(tempList) and abs(tempList[actualPosIndex] - expectedPos) <= threshold:
                            if DEBUG:
                                print("Appending ("+str(kmerIndex)+","+str(tempList[actualPosIndex])+")")
                            L1.append( (kmerIndex, tempList[actualPosIndex]) )
                    elif DEBUG:
                        print("expectedPos has already been used")
                # Make sure we meet the requirement of minimum amount of k-mers, and also check
                # that we don't consider a segment of a previously recorded LINE-1.
                if len(L1) >= minKmers: # Found one!!!
                    L1sPerChrmCount[chrm] += 1
                    if CSVOutput:
                        print("\t".join([str(x) for x in ([chrm]+L1)]))
                    else:
                        print_L1(L1, chrm, distList, CSVdir, sequence_type)
                    # Mark positions as used
                    for kmerIndex,kmerPos in L1:
                        if DEBUG: print("usedKmerPos["+str(kmerIndex)+"].add("+str(kmerPos)+")")
                        usedKmerPos[kmerIndex].add(kmerPos)

    return L1sPerChrmCount


def load_avg_lengths(CSVdir):
    """Load average FLI-L1 component lengths into L1Base2 module for GFF3 output.
    
    If an error is detected, an error string is returned.  Otherwise, nothing (None) is returned.
    """
    # We require a file ending with 'fli-l1.csv' within the CSV directory.
    csvFileList = [file for file in os.listdir(CSVdir) if file.endswith('fli-l1.csv')]
    if len(csvFileList) == 1:
        csvFilePath = CSVdir + csvFileList[0]

        if os.path.getsize(csvFilePath) > 0:
            L1CSV = l1base2.load_CSV_file(csvFilePath)
            #avg_lengths = l1base2.AvgLengths
            return True
        else:
            return "The file " + csvFilePath + " is empty!"
    elif len(csvFileList) == 0:
        return "The directory " + CSVdir + " must contain one CSV file ending with 'fli-l1.csv'!"
    else:
        return "There is more than one CSV file ending with 'fli-l1.csv'!"


def print_L1(L1, chrm, distList, CSVdir, sequence_type):
    """Print LINE-1 annotation in GFF3 format.

    Parameters:
        L1       - Pattern of hits used to detect the LINE-1
        chrm     - Chromosome in which the LINE-1 was detected
        distList - Positions of probes (to a fixed reference point)
    """
    print(chrm + "\t", end="")                # Column 1: "seqid"
    print("L1PD\t", end="")                   # Column 2: "source"
    print("mobile_genetic_element\t", end="") # Column 3: "type"

    avg_lengths = l1base2.AvgLengths if sequence_type else dbRIP.AvgLengths


    firstKmerIndex = L1[0][0]
    firstKmerPos   = L1[0][1]
    forwardStrand = False if firstKmerPos < 0 else True
    ORF1Start = int(firstKmerPos) - distList[firstKmerIndex] # Pos - distance from start of ORF1
    AVG5UTRLENGTH = avg_lengths["5UTR"]
    AVGINTERLENGTH = avg_lengths["Inter"]
    AVGORF1LENGTH = avg_lengths["ORF1"]
    AVGORF2LENGTH = avg_lengths["ORF2"]
    AVG3UTRLENGTH = avg_lengths["3UTR"]
    L1Start = abs(ORF1Start - AVG5UTRLENGTH)
    L1End = abs(ORF1Start + AVGORF1LENGTH + AVGINTERLENGTH + AVGORF2LENGTH + AVG3UTRLENGTH)

    # Start and End need to be in ascending order
    if forwardStrand:
        print(str(L1Start) + "\t", end="")# Column 4: "start"
        print(str(L1End) + "\t", end="")  # Column 5: "end"
    else:
        print(str(L1End) + "\t", end="")  # Column 4: "start"
        print(str(L1Start) + "\t", end="")# Column 5: "end"
    print(".\t", end="")              # Column 6: "score"

    print("+\t" if forwardStrand else "-\t", end="") # Column 7: "strand"
    print(".\t", end="")              # Column 8: "phase"
    print("Name=LINE1")               # Column 9: "attributes"


def main(SAMFile, fastaFile, minKmers, CSVdir, sequence_type):
    """Execute the L1PD algorithm and print the L1s that were found.

    Parameters:
        SAMFile   - File with the alignments of the probes to the subject genome
        fastaFile - FASTA file containing the probes
        minKmers  - Minimum amount of probes in pattern to consider it a LINE-1
        CSVdir    - Directory with L1Base2 metadata in CSV format
    """
    # Before we proceed any further, make sure we can access the FLI-L1 metadata.
    if sequence_type:
        result = load_avg_lengths(CSVdir)
        if result is not None and result != True:
            sys.exit(result) # Exit with error message

    # Load contents of SAM file into a dictionary where k-mer name is the key
    kmerPosDict = dict()
    # Sometimes SAM files generated by mrfast contain "invalid"
    # UTF-8 characters (not sure why), so use a different encoding.
    with open(SAMFile, encoding="ISO-8859-1") as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            pieces = line.split("\t")
            kmerName = pieces[0]
            chrm = pieces[2] # chromosome
            if len(chrm) > 5: # Only process chri; i = 1, 2, ..., 22, X, Y
                continue
            # If the FLAG field has 5th bit (16) on, seq is reverse complemented.
            # In that case we store the position as negative to indicate so.
            kmerPos = int(pieces[3]) if not int(pieces[1])&16 else -int(pieces[3])
            if chrm not in kmerPosDict:
                kmerPosDict[chrm] = dict()
            kmerPosDict[chrm][kmerName] = kmerPosDict[chrm].get(kmerName,list()) + [kmerPos]
    
    if len(kmerPosDict) == 0: # No k-mer hits?
        exit()

    # Open FASTA file to read distances from k-mers to beginning of ORF1
    distances = list()
    kmerNames = list()
    with open(fastaFile) as fh:
        for line in fh:
            if line.startswith(">"):
                pieces = line[1:].rstrip().split() # Remove beginning '>' and trailing newline first
                # First thing in description line should be the k-mer name
                kmerNames.append(pieces[0]) # Preserve order of k-mers in FASTA file
                # Last thing in description line should be the distance to start of ORF1
                distances.append(int(pieces[-1]))
    if DEBUG:
        print("kmerNames =", kmerNames)
        print("distances =", distances)

    # If minimum number of k-mers is not specified, or if the minimum specified
    # is larger than amount of k-mers, then search for all of them.
    if minKmers == 0 or abs(minKmers) >= len(kmerNames):
        minKmers = len(kmerNames)
    elif minKmers < 0: # User specified how many can be missed
        minKmers = len(kmerNames) + minKmers # Note that abs(minKmers) < len(kmerNames)

    # Sort the positions of "hits" for each k-mer.
    for chrm in kmerPosDict:
        for kmerName in kmerPosDict[chrm]:
            kmerPosDict[chrm][kmerName].sort()
            if DEBUG: print("kmerPosDict["+chrm+"]["+kmerName+"] =", kmerPosDict[chrm][kmerName])    

    L1sPerChrmCount = find_L1s(kmerPosDict, kmerNames, distances, minKmers, CSVdir, sequence_type)

    # We used to draw a histogram of the L1s found in each chromosome of the
    # user's genome vs. L1s found in each chromosome of the reference genome.
    # But, reference genome data was hard-coded for humans.
    # draw_histogram(L1sPerChrmCount)


def draw_histogram(countDict):
    """Draw a histogram (bar graph) of amount of L1s found in each chromosome.
    
    L1s found in the genome provided by the user are plotted against the
    L1s found in the reference genome.
    X axis shows the chromosomes and Y axis shows the L1 count per chromosome.
    NOTE: This script is not currently being used, but is left here for
          reference, in case future users wish to implement something similar.
    """
    # The keys have the format 'chr1', 'chr2', etc., so we remove the 'chr'
    # For efficiency, we hard-code the histogram for
    # GRCh38 with decoy (from 1000 Genome Project)
    # These counts were generated by running this script on GRCh38DH and
    # counting the amount of L1s detected for each chromosome per the GFF3.
    refHist = {'1': 721, '2': 777, '3': 800, '4': 886, '5': 779, '6': 661, '7': 513, '8': 592, '9': 401, '10': 416, '11': 547, '12': 477, '13': 250, '14': 313, '15': 239, '16': 144, '17': 89, '18': 203, '19': 88, '20': 130, '21': 59, '22': 59, 'X': 946, 'Y': 225}
    countDict.pop('chrM', None) # Remove mitochondria count if present
    xVals = [ s[3:] for s in countDict.keys() ] # Remove 'chr' from each key
    yVals = list(countDict.values())
    yValsRef = list(refHist.values())
    plt.title('Histogram of LINE-1s per chromosome')
    plt.xlabel('Chromosome')
    plt.ylabel('LINE-1 Count')
    w = 0.4
    bar1 = np.arange(len(refHist))
    bar2 = [i+w for i in bar1]
    plt.bar(bar1, yVals, w, label="Target")
    plt.bar(bar2, yValsRef, w, label="GRCh38DH")
    # We move the tick positions and use chrm names
    plt.xticks(bar1+w/2, xVals)
    plt.legend()
    plt.savefig('L1PD_hist.png')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find LINE-1s from k-mer alignments using distance from k-mers to beginning of LINE-1 consensus")
    parser.add_argument("SAM", help="Full path to SAM file with alignments")
    parser.add_argument("FASTA", help="Full path to FASTA file with ORF k-mers")
    parser.add_argument("-t", "--threshold", type=int, default=700, help="Maximum allowed difference between expected distance between k-mers and actual distance")
    parser.add_argument("-m", "--minkmers", type=int, default=2, help="Minimum amount of matched k-mers to signal a LINE-1 was found")
    parser.add_argument("--debug", help="Show debugging messages", action="store_true")
    parser.add_argument("--data_dir", help="Full path to directory with L1Base2 metadata CSV files")
    parser.add_argument("-q", type=int, default=1, help="Type of sequence: LINE or Other")
    # We need output in CSV format if determining precision & recall.
    # This will be removed in future versions, once the p&r scripts are able
    # to use GFF3 output for this purpose.
    # For now, hide the option since it's meant for internal use.
    parser.add_argument("--csvoutput", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    threshold = args.threshold
    DEBUG = args.debug
    CSVOutput = args.csvoutput
    main(args.SAM, args.FASTA, args.minkmers, args.data_dir, args.q)
