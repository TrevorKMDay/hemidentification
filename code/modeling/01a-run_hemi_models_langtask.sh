#!/bin/bash

mkdir -p results/langtask/

m=inputs/hemiconnectome_language.pickle
data=langtask

for j in A B C D E ; do

    for k in sv nn lda ; do

        echo "${j} ${k} 1"

        fname=method-${k}_outcome-hemi_test-${j}_hemi-all_data-${data}

        if [ ! -e results/${data}/${fname}.csv ] ; then

            python run_model.py                                 \
                --output_name   results/${data}/${fname}.csv    \
                class ${m} ${j} ${k}

        else

            echo "Skipping ${j} ${k} ${m} - output exists"

        fi

        fname2=method-${k}_outcome-hemi_test-${j}_hemi-all_data-${data}_subset-200

        if [ ! -e results/${data}/${fname2}.csv ] ; then

            python run_model.py                                 \
                --output_name   results/${data}/${fname2}.csv    \
                --subset        200                             \
                class ${m} ${j} ${k}

        else

            echo "Skipping ${j} ${k} ${m} - output exists"

        fi

    done

done

for i in lefty righty nonrighty righty2; do

    for j in A B C D E ; do

        for k in sv nn lda ; do

            echo "${j} ${k} 1"

            fname=method-${k}_outcome-hemi_test-${j}_hemi-all_hands-${i}_data-${data}

            if [ ! -e results/${data}/${fname}.csv ] ; then

                python run_model.py                                 \
                    --hands         ${i}                            \
                    --output_name   results/${data}/${fname}.csv    \
                    class ${m} ${j} ${k}

            else

                echo "Skipping ${i} ${j} ${k} ${m} - output exists"

            fi

            if [[ ${i} == "righty" || ${i} == "righty2" ]] ; then

                fname2=method-${k}_outcome-hemi_test-${j}_hemi-all_hands-${i}_data-${data}_subset-200

                if [ ! -e results/${data}/${fname2}.csv ] ; then

                    python run_model.py                                 \
                        --output_name   results/${data}/${fname2}.csv    \
                        --subset        200                             \
                        class ${m} ${j} ${k}

                else

                    echo "Skipping ${j} ${k} ${m} - output exists"

                fi

            fi

        done

    done

done