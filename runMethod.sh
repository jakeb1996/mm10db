#!/bin/bash

#array=( acetylcholin adenosine adrenalin calcineurin camk2 clock dopamine gaba_and_others glutamate gsk3_and_camk_others histamin prkc serotonine oldSet_part1 oldSet_part2 )
array=( acetylcholin dopamine sleep_wake adenosine adrenalin calcineurin camk2 clock gaba_and_others glutamate gsk3_and_camk_others histamin prkc serotonine )

for i in "${array[@]}"
do
    python target_identitification_viaC.py nb_threads_C=8 nb_threads_Bowtie=8 genes=$i
    cp ./targets/accepted_targets_ExcelFriendly_$i.tsv ~/Dropbox/Results/CRISPR/
done