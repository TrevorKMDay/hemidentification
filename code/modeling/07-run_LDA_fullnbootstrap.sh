#!/bin/bash

mkdir -p lda_full/

outcome=hemi
data=inputs/hemiconnectome.pickle

# One-class model (hemispheres) ===============================================

# python run_model.py                         \
#     --output_name   "lda_full/full.csv"     \
#     --hands         righty                  \
#     ${outcome} ${data} within lda

# python dump_scalings.py "lda_full/full.pickle" "lda_full/full_scalings.csv"

# for i in $(seq -w 0 9999) ; do

#     python run_model.py                                 \
#         --output_name   "lda_full/full_bs${i}.csv"      \
#         --bootstrap     1                               \
#         --hands         righty                          \
#         ${outcome} ${data} within lda

#     python dump_scalings.py \
#         "lda_full/full_bs${i}.pickle" "lda_full/full_bs${i}_scalings.csv"

# done

# One-class model (lefties) ===============================================

python run_model.py                                 \
    --output_name   "lda_full/full_lefties.csv"     \
    --hands         lefty                           \
    ${outcome} ${data} within lda

python dump_scalings.py \
    "lda_full/full_lefties.pickle" "lda_full/full_lefties_scalings.csv"

for i in $(seq -w 0 999) ; do

    python run_model.py                                         \
        --output_name   "lda_full/full_lefties_bs${i}.csv"      \
        --bootstrap     1                                       \
        --hands         righty                                  \
        ${outcome} ${data} within lda

    python dump_scalings.py \
        "lda_full/full_lefties_bs${i}.pickle" \
        "lda_full/full_lefties_bs${i}_scalings.csv"

done

# Four-class model ===============================================

# python run_model.py                                 \
#     --output_name   "lda_full/full_class.csv"       \
#     class ${data} within lda

python dump_scalings.py \
    "lda_full/full_class.pickle" \
    "lda_full/full_class_scalings.csv"

# for i in $(seq -w 0 9999) ; do

#     if [ -e "lda_full/full_class_bs${i}_scalings.csv" ] ; then continue ; fi

#     python run_model.py                                         \
#         --output_name   "lda_full/full_class_bs${i}.csv"        \
#         --bootstrap     1                                       \
#        class ${data} within lda

#     python dump_scalings.py                         \
#         "lda_full/full_class_bs${i}.pickle"         \
#         "lda_full/full_class_bs${i}_scalings.csv"

# done