#!/bin/bash

home="${HOME}/Library/CloudStorage/GoogleDrive-td758@georgetown.edu/My Drive/Projects/hemisphere_fingerprinting"

ptseries="${home}/data/test.ptseries.nii"
pconn="${home}/data/test.pconn.nii"

wb_command -cifti-correlation   \
    "${ptseries}"               \
    "${pconn}"                  \
    -fisher-z 

wb_command -cifti-convert       \
    -to-text                    \
        "${pconn}"              \
        "${pconn/nii/txt}"
