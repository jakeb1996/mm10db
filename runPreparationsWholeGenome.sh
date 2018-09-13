#!/bin/bash

FILES=./mm10_input/list_???.txt

for f in $FILES
do
    NAME=${f:13:8}
    python createListExons.py $NAME
    python prepareExonSequences.py $NAME
done
