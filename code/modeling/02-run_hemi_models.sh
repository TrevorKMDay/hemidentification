#!/bin/bash

mkdir -p results/{base,osampl}

# Single-handedness groups
for i in lefty righty ; do

    for j in A B C D E ; do

        for k in svm nn lda ; do

            for m in hemiconnectome hemiconnectome_oversampled connectome ; do

                # Rotate through the classes, except the original hand/class
                #   doesn't have a -1 suffix
                for n in "" "2" ; do

                    if [ ${m} == "hemiconnectome" ] ; then
                        data="base"
                        outcome="class${n}"
                    elif [ ${m} == "hemiconnectome_oversampled" ] ; then
                        data="osampl"
                        outcome="class${n}"
                    elif [ ${m} == "connectome" ] ; then
                        data="full"
                        outcome="handedness${n}"
                    fi

                    fname=method-${k}_outcome-hemi_test-${j}_hands-${i}_data-${data}

                    python run_model.py                                 \
                        --hands         ${i}                            \
                        --output_name   "results/${data}/${fname}.csv"  \
                        "${outcome}" ${m}.pickle ${j} ${k}

                done

            done

        done

    done

done
