CSV = None
L1sByChr = None

def __load_CSV_data(metaLstLines):
    """Convert a list of lines with metadata into a list of dictionaries."""
    global L1sByChr
    # First, convert each line into a list of strings.
    metaLstLsts = [ metaline.rstrip().split(',') for metaline in metaLstLines ]
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
    #print(L1sByChr)
    return metaLstDcts


def get_L1s_by_chr():
    """Return a dictionary of the L1s organized by chromosome (key)."""
    if L1sByChr is None:
        return None
    return L1sByChr.copy()


def load_CSV_file(fName):
    """Load CSV metadata from a local file into a list of dictionaries (of strings)."""
    global CSV
    with open(fName) as fhCSV:
        metaLstLines = fhCSV.readlines()
    CSV = __load_CSV_data(metaLstLines)
    
    # Calculate average (mode) lengths of the components for later use
    #calc_avg_lengths()
    
    return CSV.copy()
