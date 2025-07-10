#!/bin/bash

if [ ${#} -eq 1 ] ; then
    input_text=${1}
else
    echo "Supply input file!"
    exit 1
fi

template="${HOME}/Google Drive/My Drive/Projects/hemisphere_fingerprinting/data/glasserhcp/example_seg-Glasser.pscalar.nii"

echo "Starting conversion ..."

wb_command -cifti-convert                   \
    -from-text                              \
        "${input_text}"                     \
        "${template}"                       \
        "${input_text//txt/pscalar.nii}"

echo "Wrote ${input_text//txt/pscalar.nii}"
