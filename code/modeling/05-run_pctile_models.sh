#!/bin/bash

mkdir -p results/pctile/

inputs=$(find inputs/subsets/ -name "*.pickle" | sort)

for j in A B C D E ; do

    for k in sv nn lda ; do

        for m in ${inputs} ; do

            # Rotate through the classes, except the original hand/class
            #   doesn't have a -1 suffix

            outcome=class
            data=$(basename "${m}" .pickle | sed 's/hemiconnectome_//')

            echo "${j} ${k} ${m} ${n:=1}"

            fname=method-${k}_outcome-class_test-${j}_hands-all_data-${data}

            python run_model.py                                 \
                --output_name   "results/pctile/${fname}.csv"  \
                "${outcome}" "${m}" ${j} ${k}

        done

    done

done