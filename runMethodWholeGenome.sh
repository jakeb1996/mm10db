#!/bin/bash

FILES=./mm10_input/list_???.txt
LOG=~/Dropbox/Results/CRISPR/logTestWholeGenome.txt

rm -f $LOG

echo "Starting..." >> $LOG

for f in $FILES
do
    NAME=${f:13:8}
    date >> $LOG
    echo $NAME >> $LOG
    python target_identitification_viaC.py nb_threads_C=128 nb_threads_Bowtie=8 genes=$NAME
    cp ./targets/accepted_targets_ExcelFriendly_$NAME.tsv ~/Dropbox/Results/CRISPR/
    date >> $LOG
done

echo "Done." >> $LOG
