#!/bin/bash

mkdir -p results/langtask/

m=inputs/hemiconnectome_language.pickle
data=langtask

for j in A B C D E ; do

    for k in sv nn lda ; do

        echo "${j} ${k} 1"

        fname=method-${k}_outcome-hand_test-${j}_hemi-all_data-${data}
        fname2=method-${k}_outcome-hand2_test-${j}_hemi-all_data-${data}

        if [ ! -e results/${data}/${fname}.csv ] ; then

            python run_model.py                                 \
                --output_name   results/${data}/${fname}.csv    \
                class ${m} ${j} ${k}

        else

            echo "Skipping ${j} ${k} ${m} - output exists"

        fi

        echo "${j} ${k} 2"

        if [ ! -e results/${data}/${fname2}.csv ] ; then

            python run_model.py                                 \
                --output_name   results/${data}/${fname2}.csv    \
                class2 ${m} ${j} ${k}

        else

            echo "Skipping ${j} ${k} ${m} - output exists"

        fi

    done

done


# Single-hemisphere groups
# for i in LH RH ; do

#     for j in A B C D E ; do

#         for k in svm nn lda ; do

#             m=hemiconnectome_language.pickle
#             data=langtask

#             fname=method-${k}_outcome-hand_test-${j}_hemi-${i}_data-${data}

#             if [ ! -e results/${data}/${fname}.csv ] ; then

#                 python run_model.py                                 \
#                     --hemi          ${i}                            \
#                     --output_name   results/${data}/${fname}.csv    \
#                     class ${m} ${j} ${k}

#             else

#                 echo "Skipping ${i} ${j} ${k} ${m} - output exists"

#             fi

#         done

#     done

# done

