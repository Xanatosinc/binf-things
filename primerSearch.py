#!/usr/bin/env python
"""
primerSearch.py:  

Input: FASTA file, primer file (5'\n3')
Output: Amplicon(s) and info
"""
import re
import sys

DEFAULT_DEBUG = 0

def generateRegex(pattern):
    '''
        generateRegex(pattern) -> regexPattern from IUPAC ambiguity codes
    '''
    ambiguity = {
        'A':'A','C':'C','G':'G','T':'T',
        'U':'U','R':'[AG]','Y':'[CT]','S':'[CG]',
        'W':'[AT]','K':'[GT]','M':'[AC]',
        'B':'[CGT]','V':'[ACG]','D':'[AGT]','H':'[ACT]'
    }
    pattern = list(pattern)
    pattern = [ambiguity.get(item,item) for item in pattern]
    return ''.join(pattern)

def getFastaSeqs(ifn):
    '''
        getFastaSeqs(input file name) -> [['>seqInfo1', 'ACTG...'], ['>seqInfo2', 'CATG...']]
    '''
    seqs=[]
    seqName=""
    nucleotides=""
    with open(ifn,"r") as infile:
        for line in infile:
            if (line.startswith(">")):
                seqName = line.strip()
                if (not nucleotides==""):
                    seqs.append([oldSeqName, nucleotides])
                    nucleotides=""
                oldSeqName = seqName
            else:
                nucleotides += treatUracil(line.strip().upper())
        seqs.append([oldSeqName, nucleotides])
    return seqs

def treatUracil(seq):
    '''
        treatUracil(sequence) -> swap Us for Ts
    '''
    uracil = {'U':'T'}
    seq = list(seq)
    seq = [uracil.get(item,item) for item in seq]
    return ''.join(seq)

def getPrimers(ifn):
    '''
        getPrimers(inputFileName) -> ['ACCYG...', 'YGGRC...']
    '''
    validPattern = r'[^\.ACGTUMRWSYKVHDB]'
    with open(ifn,"r") as infile:
        lines = infile.readlines()
        threePrime = lines[0].strip().upper()
        fivePrime = lines[1].strip().upper()
        if (re.search(validPattern, threePrime) or re.search(validPattern, fivePrime)):
            print("Primer(s) contain invalid character(s). Limit to [ACGTUMRWSYKVHDB] or take it up with IUPAC.")
            sys.exit(1)
        threePrime = generateRegex(threePrime)
        fivePrime = generateRegex(fivePrime)
        sys.exit

    return [threePrime, fivePrime]

def getSeqInd(seq, pattern):
    '''
        getSeqInd(sequence, regex pattern) -> list of indeces
    '''
    ind = []                            # Initialize list of indeces
    lastind = 0                         # Previous index
    pre = re.compile(pattern)
    while (pre.search(seq, lastind)):
        lastind = pre.search(seq, lastind).start()
        ind.append(lastind)
        lastind += 1
    return list(ind)

def insult(sequence):
    '''
        insult(sequence) -> reverse complement of sequence
    '''
    complement = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
        'U': 'A', 'R': 'Y', 'Y': 'R', 'S': 'S',
        'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V',
        'V': 'B', 'D': 'H', 'H': 'D',
    }
    bases = list(sequence)
    bases = reversed([complement.get(base,base) for base in bases])
    return ''.join(bases)

def revInd(indices, factor):
    '''
        revInd(iterable indices, int factor) -> [factor-ind1, factor-ind2...]
    '''
    rev = []
    for i in range(len(indices)):
        rev.append(factor-indices[i-1])
    return rev

# Command line arguments
args = sys.argv
debugActive = DEFAULT_DEBUG
if ('-d' in args):
    debugActive = 1
    args.remove('-d')
if (len(args) != 3 or '-h' in args):
    print "Usage: primerSearch.py [-d] FASTA PRIMERS"
    print "PRIMERS is a two-line file with ambiguity-code primers."
    print "-d Debug mode: Include entire fasta sequence [U->T] in output."
    sys.exit(1)
else:
	fastaFileName = sys.argv[1] # Define FASTA input file name
	primerFileName = sys.argv[2] # Define primer input file name
	
seqPairs = getFastaSeqs(fastaFileName)          # get sequence from file
primers = getPrimers(primerFileName)            # get primers from file

for seqPair in seqPairs:
    seqName = seqPair[0]
    sequence = seqPair[1]

    startList = getSeqInd(sequence, primers[0])
    endListRev = getSeqInd(insult(sequence), primers[1]) # Low value means end of string
    endList = revInd(endListRev, len(sequence)) # Low value means beginning of string

    start = min(startList) if len(startList) else -1
    endFp = max(endList) if len(endList) else -1 # chars from beginning of sequence at which to stop
    endTp = min(endListRev) if len(endListRev) else -1 # chars from terminus of sequence at which to stop

    # Print things
    if (debugActive):
        print sequence
    if (start > -1 and endFp > -1):
        print(seqName)
        if (endFp > start):
            amplicon = sequence[start:-endTp]
            friendlyStarts = [item+1 for item in startList]
            friendlyEnds = [item+1 for item in endList]
            print(["5'", primers[0], "3'", primers[1], "starts:", friendlyStarts, "ends:", friendlyEnds, "amplicon len:", len(amplicon)])
            print(amplicon)
        else:
            print "Note: start and end out of order."
    else:
        print "Note: no primer match."
        friendlyStarts = [item+1 for item in startList]
        friendlyEnds = [item+1 for item in endList]
        print(["5'", primers[0], "3'", primers[1], "starts:", friendlyStarts, "ends:", friendlyEnds, "amplicon len:", len(amplicon)])
        
sys.exit
