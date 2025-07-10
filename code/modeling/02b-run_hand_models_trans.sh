#!/bin/bash

mkdir -p results/trans/

# Single-handedness groups
data=inputs/transconnectome.pickle

for j in A B C D E ; do

    for k in sv nn lda ; do

        fname=method-${k}_outcome-hand_test-${j}_data-trans

        python run_model.py                         \
            --output_name   results/trans/${fname}.csv    \
            handedness ${data} ${j} ${k}

    done

done

for j in A B C D E ; do

    for k in sv nn lda ; do

        fname=method-${k}_outcome-hand2_test-${j}_data-trans

        python run_model.py                         \
            --output_name   results/trans/${fname}.csv    \
            handedness2 ${data} ${j} ${k}

    done

done