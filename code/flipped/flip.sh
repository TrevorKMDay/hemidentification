#!/bin/bash

mirror=atlas-GlasserMirrored_space-fsLR_den-32k_dseg.dscalar.nii

if [ ! -e ${mirror} ] ; then

    ./rois_create_mirror.sh \
        "$(which wb_command)"                                       \
        atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii            \
        ${mirror}

fi

wb_command -cifti-label-export-table \
    atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii    \
    1                                                   \
    atlas-Glasser_labels.txt

sed -e 's/^R_/R2L_/' -e 's/^L_/L2R_/' atlas-Glasser_labels.txt > \
    atlas-GlasserMirrored_labels.txt

wb_command -cifti-label-import          \
    ${mirror}                           \
    atlas-GlasserMirrored_labels.txt    \
    ${mirror/dscalar/dlabel}