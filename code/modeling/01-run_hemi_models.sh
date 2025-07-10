#!/bin/bash

mkdir -p results/{base,osampl,full}

# for j in A B C D E ; do

#     for k in sv nn lda ; do

#         for m in hemiconnectome hemiconnectome_oversampled connectome ; do

#             # Rotate through the classes, except the original hand/class
#             #   doesn't have a -1 suffix
#             for n in "" "2" ; do

#                 if [ ${m} == "hemiconnectome" ] ; then
#                     data="base"
#                     outcome="class${n}"
#                 elif [ ${m} == "hemiconnectome_oversampled" ] ; then
#                     data="osampl"
#                     outcome="class${n}"
#                 elif [ ${m} == "connectome" ] ; then
#                     data="full"
#                     outcome="handedness${n}"
#                 fi

#                 echo "${j} ${k} ${m} ${n:=1}"

#                 fname=method-${k}_outcome-hemi_test-${j}_hands-all_data-${data}

#                 python run_model.py                                 \
#                     --output_name   "results/${data}/${fname}.csv"  \
#                     "${outcome}" inputs/${m}.pickle ${j} ${k}

#             done

#         done

#     done

# done

# # Single-handedness groups
# for i in lefty righty ; do

#     for j in A B C D E ; do

#         for k in sv nn lda ; do

#             for m in hemiconnectome hemiconnectome_oversampled ; do

#                 # Rotate through the classes, except the original hand/class
#                 #   doesn't have a -1 suffix
#                 for n in "" "2" ; do

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

#                     echo "${i} ${j} ${k} ${m} ${n}"

#                     fname=method-${k}_outcome-hemi_test-${j}_hands-${i}_data-${data}

#                     python run_model.py                                 \
#                         --hands         ${i}                            \
#                         --output_name   "results/${data}/${fname}.csv"  \
#                         "${outcome}" inputs/${m}.pickle ${j} ${k}

#                 done

#             done

#         done

#     done

# done

for i in nonrighty righty2 ; do

    for j in A B C D E ; do

        for k in sv nn lda ; do

            for m in hemiconnectome hemiconnectome_oversampled ; do

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

                    echo "${i} ${j} ${k} ${m} ${n}"

                    fname=method-${k}_outcome-hemi_test-${j}_hands-${i}_data-${data}

                    python run_model.py                                 \
                        --hands         ${i}                            \
                        --output_name   "results/${data}/${fname}.csv"  \
                        "${outcome}" inputs/${m}.pickle ${j} ${k}

                        # break 5

                done

            done

        done

    done

done
