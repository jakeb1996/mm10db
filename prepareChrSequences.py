
##################################
##                              ##
##    prepareChrSequences.py    ##
##                              ##
##################################

#
# Purpose: preparing all the chromosome sequences
#

#
# Inputs:
#           - one .fa file for each chromosome
#

#
# Output:
#           - one file with all the sequences, one chromosome per line
#


from time import localtime, strftime, sleep

dir = "./mm10_input/chr_sequences/"

chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "M", "X", "Y"]

outFile = open(dir+"all_sequences.txt",'w')

for c in chromosomes:
    chr = "chr"+c+".fa"
    print strftime("%H:%M:%S", localtime())+": "+chr[:-3]
    inFile = open(dir+chr,'r')
    lines = inFile.readlines()
    for i in range(1,len(lines)):
        outFile.write(lines[i].rstrip())
    outFile.write("\n")

outFile.close()

print strftime("%H:%M:%S", localtime())+": Done."



#######################
##                   ##
##    End of File    ##
##                   ##
#######################
