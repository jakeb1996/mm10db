#!/bin/bash

array=( acetylcholin adenosine adrenalin calcineurin camk2 clock dopamine gaba_and_others glutamate gsk3_and_camk_others histamin prkc serotonine oldSet_part1 oldSet_part2 sleep_wake )

for i in "${array[@]}"
do
    python createListExons.py $i
    python prepareExonSequences.py $i
done
