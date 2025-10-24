#!/bin/bash

mkdir -p results/n24/

subset=${1}

# Single-handedness groups
i=righty
m=hemiconnectome
data=base
# subset=24

for j in A B C D E ; do

    for k in sv nn lda ; do

        # Rotate through the classes, except the original hand/class
        #   doesn't have a -1 suffix
        for n in "" "2" ; do

            outcome="class${n}"

            echo "${i} ${j} ${k} ${m} hand${n}"

            fname=method-${k}_outcome-hemi_test-${j}_hands-${i}_data-${data}_subset-${subset}

            python run_model.py                                 \
                --hands         ${i}                            \
                --output_name   "results/n24/${fname}.csv"      \
                --subset        "${subset}"                     \
                "${outcome}" inputs/${m}.pickle ${j} ${k}

            # break 3

        done

    done

done

i=righty2

for j in A B C D E ; do

    for k in sv nn lda ; do

        # Rotate through the classes, except the original hand/class
        #   doesn't have a -1 suffix
        for n in "" "2" ; do

            outcome="class${n}"

            echo "${i} ${j} ${k} ${m} ${n}"

            fname=method-${k}_outcome-hemi_test-${j}_hands-${i}_data-${data}_subset-${subset}

            python run_model.py                                 \
                --hands         ${i}                            \
                --output_name   "results/n24/${fname}.csv"      \
                --subset        "${subset}"                     \
                "${outcome}" inputs/${m}.pickle ${j} ${k}

            # break 3

        done

    done

done