#!/bin/bash

mkdir -p lda_full/

outcome=hemi
data=inputs/hemiconnectome.pickle

# One-class model (hemispheres) ===============================================

# This runs the hemisphere classification model on all righties, and the
#   associated bootstrap

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

# This runs the hemisphere classification model on all lefties, and the
#   associated bootstrap

# python run_model.py                                 \
#     --output_name   "lda_full/full_lefties.csv"     \
#     --hands         lefty                           \
#     ${outcome} ${data} within lda

# python dump_scalings.py \
#     "lda_full/full_lefties.pickle" "lda_full/full_lefties_scalings.csv"

# for i in $(seq -w 0 999) ; do

#     python run_model.py                                         \
#         --output_name   "lda_full/full_lefties_bs${i}.csv"      \
#         --bootstrap     1                                       \
#         --hands         righty                                  \
#         ${outcome} ${data} within lda

#     python dump_scalings.py \
#         "lda_full/full_lefties_bs${i}.pickle" \
#         "lda_full/full_lefties_bs${i}_scalings.csv"

# done

# Four-class model ===============================================

# python run_model.py                                 \
#     --output_name   "lda_full/full_class.csv"       \
#     class ${data} within lda

# python dump_scalings.py \
#     "lda_full/full_class.pickle" \
#     "lda_full/full_class_scalings.csv"

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

# Four-class model (run 1) ==========

python run_model.py                                      \
    --output_name   "lda_full/full_class_run1.csv"       \
    class inputs/half_hemiconnectome.pickle within lda

# args: model data out
python fit_transform.py \
    lda_full/full_class_run1.pickle     \
    inputs/half_hemiconnectome.pickle   \
    lda_full/run1_xform.csv

# Four-class model (Z-scored) ==========

python run_model.py                                     \
    --output_name   "lda_full/full_class_zscored.csv"   \
    class inputs/hemiconnectome_Z.pickle within lda

python dump_scalings.py \
    "lda_full/full_class_zscored.pickle" \
    "lda_full/full_class_zscored_scalings.csv"
