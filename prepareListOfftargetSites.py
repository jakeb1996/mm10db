
####################################
##                                ##
##  prepareListOfftargetSites.py  ##
##                                ##
####################################

#
# Purpose: identify all offtarget sites in the whole genome
#

#
# Inputs:
#           - one file with all the chromosome sequences, one chromosome per line
#

#
# Output:
#           - one file with all the sites
#


from time import localtime, strftime
import re
import string


# Defining the patterns used to detect sequences
pattern_forward_offsite = r"(?=([ACG][ACGT]{19}[ACGT][AG]G))"
pattern_reverse_offsite = r"(?=(C[CT][ACGT][ACGT]{19}[TGC]))"


dir_seq = "./mm10_input/chr_sequences/"

inputFile_chr = "all_sequences.txt"
outputFile = "offtargetSites.txt"

chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "M", "X", "Y"]


# Function that returns the reverse-complement of a given sequence
def rc(dna):
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = dna.translate(complements)[::-1]
    return rcseq





inFile = open(dir_seq+inputFile_chr,'r')
outFile = open(dir_seq+outputFile,'w')

lineNumber = 0

# For every chromosome
for line_chr in inFile:

    print strftime("%H:%M:%S", localtime())+":\tChromomose "+chromosomes[lineNumber]
    
    # we parse the line and look for forward sequences
    print strftime("%H:%M:%S", localtime())+":\t\tForward-parsing the chromosome."
    match_chr = re.findall(pattern_forward_offsite,line_chr)
    print strftime("%H:%M:%S", localtime())+":\t\tProcessing and saving the parsed sequences."
    if match_chr:
        print "\t\t\tWe detected "+str(len(match_chr))+" possible off-target sites."
        # we save each sequence
        for i in range(0,len(match_chr)):
            outFile.write(match_chr[i][0:20]+"\n")
    else:
        print "\t\t\tWe did not detect any possible off-target sites."

    # we parse the line and look for reverse sequences
    print strftime("%H:%M:%S", localtime())+":\t\tReverse-parsing the chromosome."
    match_chr = re.findall(pattern_reverse_offsite,line_chr)
    print strftime("%H:%M:%S", localtime())+":\t\tProcessing the parsed sequences."
    if match_chr:
        print "\t\t\tWe detected "+str(len(match_chr))+" possible off-target sites."
        # we save each reverse-complement sequence
        for i in range(0,len(match_chr)):
            # we reverse-complement the sequence and count the number of mismatches between the rc and the target
            outFile.write(rc(match_chr[i])[0:20]+"\n")

    lineNumber += 1

inFile.close()
outFile.close()


print "\n"+strftime("%H:%M:%S", localtime())+":\tDone."
