#!/bin/bash

mkdir -p results/ehi/
outcome="EHI"

for j in A B C D E ; do

    for k in sv nn ; do

        for m in hemiconnectome hemiconnectome_oversampled connectome ; do

            # Rotate through the classes, except the original hand/class
            #   doesn't have a -1 suffix

            if [ ${m} == "hemiconnectome" ] ; then
                data="base"
            elif [ ${m} == "hemiconnectome_oversampled" ] ; then
                data="osampl"
            elif [ ${m} == "connectome" ] ; then
                data="full"
            fi

            fname="method-${k}_outcome-ehi_test-${j}_hemi-all_data-${data}"

            if [ ! -e "results/ehi/${fname}.csv" ] ; then

                python run_model.py                                 \
                    --output_name   "results/ehi/${fname}.csv"  \
                    "${outcome}" ${m}.pickle ${j} ${k}

                # break 3

            else

                echo "Skipping ${j} ${k} ${m} - output exists"

            fi

        done

    done

done