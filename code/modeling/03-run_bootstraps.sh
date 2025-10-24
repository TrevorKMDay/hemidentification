#!/bin/bash

input=inputs/hemiconnectome.pickle
group=A

for i in sv lda nn ; do

    python run_model.py                 \
        --output_name   "testing_${i}.csv"   \
        class ${input} ${group} ${i}

done

for i in sv lda nn ; do

    for j in lefty righty nonrighty righty2 ; do

        python run_model.py                 \
            --hands ${j}                        \
            --output_name   "testing_${i}_${j}.csv"   \
            class ${input} ${group} ${i}

    done

done