#!/bin/bash

if [ "${CONDA_DEFAULT_ENV}" != "hemiconnectome" ] ; then
    echo "Activate correct environment:"
    echo "  conda activate hemiconnectome"
    exit 1
fi

mkdir -p results/pctile

for j in A B C D E ; do

    for k in sv nn lda ; do

        for pctile in 0.5 0.6 0.7 0.8 0.9 0.95 0.975 0.99 ; do

            m="hemiconnectome_pctile${pctile}.pickle"
            data="p${pctile}"
            outcome1="hemi"
            outcome2="class"

            fname1="method-${k}_outcome-${outcome1}_test-${j}_data-${data}"
            fname2="method-${k}_outcome-${outcome2}_test-${j}_data-${data}"

            if [ ! -e "results/pctile/${fname1}.csv" ] ; then

                python run_model.py                                 \
                    --output_name   "results/pctile/${fname1}.csv"  \
                    "${outcome1}" ${m} ${j} ${k}

                # break 3

            else

                echo "Skipping ${j} ${k} ${m} - output exists"

            fi

            if [ ! -e "results/pctile/${fname2}.csv" ] ; then

                python run_model.py                                 \
                    --output_name   "results/pctile/${fname2}.csv"  \
                    "${outcome2}" ${m} ${j} ${k}

                # break 3

            else

                echo "Skipping ${j} ${k} ${m} - output exists"

            fi

        done

    done

done
