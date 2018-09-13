
##################################
##                              ##
##   prepareExonSequences.py    ##
##                              ##
##################################

#
# Purpose: preparing the sequences for all exons of interest
#

#
# Inputs:
#           - one file with all the chromosome sequences, one chromosome per line
#           - the list of exons of interest
#

#
# Output:
#           - one file with all the exon sequences, one exon per line
#           - note: we include 17bp before and after the exon
#


from time import localtime, strftime, sleep
import re
import ast
import sys
import os


dir_seq = "./mm10_input/chr_sequences/"
dir_list = "./mm10_input/"
chr_ = "all_sequences.txt"
list_ = "exonList.txt"
exons_ = "exon_sequences.txt"


if len(sys.argv) > 2:
    print "Usage: python prepareExonSequences.py [ID for list of genes]"
    quit ()

if len(sys.argv) == 1:
    print "Extracting exon sequences for the whole genome"

else:
    list_ = list_[:-4]+"_"+sys.argv[1]+".txt"
    print "Extracting sequence for exons listed in "+list_
    if os.path.isfile(dir_list+list_)==False:
        print "Error: file "+list_+" does not exist."
        quit()
    exons_ = exons_[:-4]+"_"+sys.argv[1]+".txt"


padding = 17 # digestion site must be on exon, but the whole target does not have to

chr_offset = {"1": 0, "2": 1, "3": 2, "4": 3, "5": 4, "6": 5, "7": 6, "8": 7, "9": 8, "10": 9, "11": 10, "12": 11, "13": 12, "14": 13, "15": 14, "16": 15, "17": 16, "18": 17, "19": 18, "M": 19, "X": 20, "Y": 21}

pattern = r".+\tchr(.+)\t[\+-]\t(\d+)\t(\d+)"

# First, we load all the chromosome sequences
print strftime("%H:%M:%S", localtime())+": Loading chromosome sequences from "+chr_
inFile = open(dir_seq+chr_,'r')
chrSeq = inFile.readlines()
inFile.close()


# Then, we load all the list of exons
print strftime("%H:%M:%S", localtime())+": Loading the exon list from "+list_
inFile = open(dir_list+list_,'r')
exonList = inFile.readlines()
inFile.close()


# Finally, we can extract (and save) the sequence for each exon (with padding on each side)
print strftime("%H:%M:%S", localtime())+": Extracting and saving the exon sequences"
outFile = open(dir_seq+exons_,'w')
for i in range(0,len(exonList)):
    match = re.search(pattern,exonList[i])
    if match:
        chr = chr_offset[match.group(1)]
        start = ast.literal_eval(match.group(2))-1
        end = ast.literal_eval(match.group(3))
        outFile.write(chrSeq[chr][start-padding:end+padding].upper()+"\n")
    else:
        print "Problem? "+exonList[i].rstrip()


print strftime("%H:%M:%S", localtime())+": Done."



#######################
##                   ##
##    End of File    ##
##                   ##
#######################
