#!/bin/bash

mkdir -p results/strong_lefties/
outcome="hemi"

for j in $(seq -w 01 29) ; do

    for k in sv nn lda ; do

        for m in strong_lefties ; do

            # Rotate through the classes, except the original hand/class
            #   doesn't have a -1 suffix

            fname="method-${k}_outcome-ehi_test-${j}_hemi-all_data-stronglefties"

            if [ ! -e "results/strong_lefties/${fname}.csv" ] ; then

                python run_model.py                                         \
                    --output_name   "results/strong_lefties/${fname}.csv"   \
                    "${outcome}" inputs/${m}.pickle "${j}" ${k}

                # break 3

            else

                echo "Skipping ${j} ${k} ${m} - output exists"

            fi

        done

    done

done