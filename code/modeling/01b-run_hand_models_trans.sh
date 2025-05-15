#!/bin/bash

mkdir -p results/

# Single-handedness groups
data=transconnectome.pickle

for j in A B C D E ; do

    for k in svm nn lda ; do

        fname=method-${k}_outcome-hand_test-${j}_data-trans

        python run_model.py                         \
            --output_name   results/${fname}.csv    \
            handedness ${data} ${j} ${k}

    done

done
