#!/bin/bash

if [ "${CONDA_DEFAULT_ENV}" != "hemiconnectome" ] ; then
    echo "Activate correct environment:"
    echo "  conda activate hemiconnectome"
    exit 1
fi

mkdir -p results/{base,osampl,full}

for j in A B C D E ; do

    for k in sv nn lda ; do

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

                fname="method-${k}_outcome-hand${n}_test-${j}_hemi-all_data-${data}"

                if [ ! -e "results/${data}/${fname}.csv" ] ; then

                    python run_model.py                                 \
                        --output_name   "results/${data}/${fname}.csv"  \
                        "${outcome}" inputs/${m}.pickle ${j} ${k}

                    # break 3

                else

                    echo "Skipping ${j} ${k} ${m} - output exists"

                fi

            done

        done

    done

done

# Single-hemisphere groups
# for i in LH RH ; do

#     for j in A B C D E ; do

#         for k in svm nn lda ; do

#             for m in hemiconnectome hemiconnectome_oversampled connectome ; do

#                 # Rotate through the classes, except the original hand/class
#                 #   doesn't have a -1 suffix
#                 for n in "" "2" ; do

#                     echo "${k} ${m} ${i} ${j} class${n}"

#                     if [ ${m} == "hemiconnectome" ] ; then
#                         data="base"
#                         outcome="class${n}"
#                     elif [ ${m} == "hemiconnectome_oversampled" ] ; then
#                         data="osampl"
#                         outcome="class${n}"
#                     elif [ ${m} == "connectome" ] ; then
#                         data="full"
#                         outcome="handedness${n}"
#                     fi

#                     fname=method-${k}_outcome-hand_test-${j}_hemi-${i}_data-${data}

#                     if [ ! -e "results/${data}/${fname}.csv" ] ; then

#                         python run_model.py                                 \
#                             --hemi          ${i}                            \
#                             --output_name   "results/${data}/${fname}.csv"  \
#                             "${outcome}" ${m}.pickle ${j} ${k}

#                     else

#                         echo "Skipping ${i} ${j} ${k} ${m} - output exists"

#                     fi

#                 done

#             done

#         done

#     done

# done

