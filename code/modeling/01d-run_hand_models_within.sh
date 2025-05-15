#!/bin/bash

if [ "${CONDA_DEFAULT_ENV}" != "hemiconnectome" ] ; then
    echo "Activate correct environment:"
    echo "  conda activate hemiconnectome"
    exit 1
fi

mkdir -p results/within

for k in sv nn lda ; do

    m="hemiconnectome"
    data="base"
    outcome="class"

    fname="method-${k}_outcome-${outcome}_test-within_hemi-all_data-${data}"

    if [ ! -e "results/${data}/${fname}.csv" ] ; then

        python run_model.py                                 \
            --output_name   "results/within/${fname}.csv"  \
            "${outcome}" ${m}.pickle within ${k}

        # break 3

    else

        echo "Skipping within ${k} ${m} - output exists"

    fi

done

