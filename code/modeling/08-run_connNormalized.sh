#!/bin/bash

if [[ ${CONDA_DEFAULT_ENV} != "hemiconnectome" ]] ; then

    echo "Set conda env to 'hemiconnectome'. Exiting."
    exit 1
fi

# This runs the hemiconnectome that has been normalized per connection
# (not across hemiconnectomes)

mkdir -p results/connNormal/

m=hemiconnectome_connNormal
data=connNormal

for j in A B C D E ; do

    for k in sv nn lda ; do

        echo "${j} ${k} ${m}"

        # One-way models
        outcome=hemi

        fnameR=method-${k}_outcome-hemi_test-${j}_hands-righty_data-${data}
        fnameL=method-${k}_outcome-hemi_test-${j}_hands-lefty_data-${data}

        python run_model.py                                 \
            --hands righty                                  \
            --output_name   "results/${data}/${fnameR}.csv" \
            "${outcome}" inputs/${m}.pickle ${j} ${k}

        python run_model.py                                 \
            --hands lefty                                   \
            --output_name   "results/${data}/${fnameL}.csv" \
            "${outcome}" inputs/${m}.pickle ${j} ${k}

        # Two-way models
        outcome=class


        fname=method-${k}_outcome-class_test-${j}_hands-all_data-${data}

        python run_model.py                                 \
            --output_name   "results/${data}/${fname}.csv"  \
            "${outcome}" inputs/${m}.pickle ${j} ${k}

    done

    # Uncomment this to run all vars for only test group A
    # break

done
