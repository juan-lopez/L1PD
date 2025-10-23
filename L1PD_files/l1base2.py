"""Module for manipulating L1Base 2.

A module for L1Base2 is created with methods to import the data and
retrieve it in a list, either by ID or by chromosome.
"""

import statistics

from io import StringIO
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from Bio import SeqIO

L1s = None # List of Bio.SeqRecord; index will be one less than the L1 ID
CSV = None # List of dictionaries of strings (e.g. CSV[0]["ORF1 Start"])
L1sByChr = None # Dictionary of L1 IDs; key will be the chromosome (e.g. 21 or X)
AvgLengths = {} # Dictionary of average (mode) lengths of components and L1s


def __retrieve_web_data(frmt):
	"""Retrieve L1Base2 file from the web.

	Parameter:
		frmt - Format in which the file should be formatted
		       1 = Fasta, 2 = GenBank, 3 = CSV
	Returns:
		Retrieved data as a string
	"""
	url = 'http://l1base.charite.de/exportall.php?DBN=hsflil1_8438&START=0&NUM=50&SORT=0'
	# format: 1 = Fasta, 2 = GenBank, 3 = CSV
	# dest: 1 = File, 2 = Window
	data = {'format':frmt, 'dest':1}
	request = Request(url, urlencode(data).encode())
	response = urlopen(request).read().decode()
	return response.strip()


def load_L1s_file(fName):
	"""Load local FASTA file and return its contents as a list of SeqRecords."""
	global L1s
	L1s = list(SeqIO.parse(fName, "fasta"))
	return L1s.copy()


def load_L1s_web():
	"""Load FASTA file from L1Base2 website and return its contents as a list of SeqRecords."""
	global L1s
	L1s = list(SeqIO.parse(StringIO(__retrieve_web_data(1)),'fasta'))
	return L1s.copy()


def load_CSV_web():
	"""Load CSV metadata from L1Base2 website into a list of dictionaries (of strings)."""
	global CSV
	metaLstLines = __retrieve_web_data(3).split("\n")
	CSV = __load_CSV_data(metaLstLines)
	return CSV.copy()


def load_CSV_file(fName):
	"""Load CSV metadata from a local file into a list of dictionaries (of strings)."""
	global CSV
	with open(fName) as fhCSV:
		metaLstLines = fhCSV.readlines()
	CSV = __load_CSV_data(metaLstLines)
	
	# Calculate average (mode) lengths of the components for later use
	calc_avg_lengths()
	
	return CSV.copy()


def __load_CSV_data(metaLstLines):
	"""Convert a list of lines with metadata into a list of dictionaries."""
	global L1sByChr
	# First, convert each line into a list of strings.  Columns are enclosed
	# in quotes and separated by a comma, except for the first and last
	# columns, so we just remove the beginning and ending quotes with [1:-1].
	metaLstLsts = [ metaline[1:-1].split('","') for metaline in metaLstLines ]
	# We'll use metaLstLsts (list of lists of strings) to create metaLstDcts
	# (list of dictionaries of strings).  The first row of the CSV contains the
	# column names, which will become the dictionary keys.
	metaLstDcts = list()
	for i in range(1,len(metaLstLsts)):
		metaDict = dict()
		for j in range(len(metaLstLsts[0])):
			metaDict[metaLstLsts[0][j]] = metaLstLsts[i][j]
		metaLstDcts.append(metaDict)

	# Since we're loading CSV data, take the opportunity to create dict of
	# L1s ordered by chromosome.  That way it's ready for when it's needed
	L1sByChr = dict()
	for i in range(len(metaLstDcts)):
		L1sByChr[metaLstDcts[i]["Chr"]] = L1sByChr.get(metaLstDcts[i]["Chr"],list()) + [i]

	return metaLstDcts


def get_L1s():
	"""Return a copy of the L1s (list of Bio.SeqRecord)"""
	return L1s.copy()


def get_L1s_by_chr():
	"""Return a dictionary of the L1s organized by chromosome (key)."""
	if L1sByChr is None:
		return None
	return L1sByChr.copy()


def get_CSV():
	"""Return a copy of the CSV metadata (list of dictionaries of strings)."""
	return CSV.copy()


def calc_avg_lengths():
	"""Calculate the average length for each L1 component and store results in dictionary."""
	global AvgLengths

	AvgLengths["5UTR"] = calc_avg_length("ORF1 Start", "", -1)
	AvgLengths["ORF1"] = calc_avg_length("ORF1 End", "ORF1 Start", 1)
	AvgLengths["Inter"] = calc_avg_length("ORF2 Start", "ORF1 End", -1)
	AvgLengths["ORF2"] = calc_avg_length("ORF2 End", "ORF2 Start", 1)
	AvgLengths["3UTR"] = calc_avg_length("ORF1 Start", "ORF2 End", 1)
	AvgLengths["L1"] = calc_avg_length("End", "Start", 1)


def calc_avg_length(keyEnd, keyStart, offset):
    """Calculate the average length for an L1 component using the CSV metadata.

    This is a helper function for calc_avg_lengths.

    Parameters:
        keyEnd   - Dictionary key to use for end position
        keyStart - Dictionary key to use for start position
        offset   - Additional value to add (usually either -1, 0, or +1)
    """
    global CSV

    values = []  # List of lengths used to determine mode
    if keyStart == "ORF2 End":  # Special case for 3'UTR
        for row in CSV:
            if row["ORF2 End"].isdigit():
                values.append(abs(abs(int(row["End"]) - int(row["Start"])) - int(row["ORF2 End"])) + 1)
    elif keyStart == "":  # Disregard this column
        for row in CSV:
            values.append(int(row[keyEnd]) + int(offset))
    else:
        for row in CSV:
            if row[keyEnd].isdigit():
                values.append(abs(int(row[keyEnd]) - int(row[keyStart])) + int(offset))

    try:
        mode_value = statistics.mode(values)
    except statistics.StatisticsError:
        # Handle multiple modes manually
        modes = [value for value in set(values) if values.count(value) == max([values.count(v) for v in set(values)])]
        mode_value = max(modes)

    return mode_value
