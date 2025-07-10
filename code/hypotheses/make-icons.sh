#!/bin/bash

data=../data/

same_hemi=${data}/brain.png
diff_hemi=${data}/arrows.png

same_hand=${data}/pray.png
diff_hand=${data}/open_hands.png

overlay() {

    magick \
        \( -page +0+0 "${1}" \)   \
        \( -page +60+60 "${2}" \) \
        -layers mosaic             \
        -fuzz 5%                   \
        -transparent white         \
        "${3}.png"

}

mkdir -p icons/

overlay ${same_hemi} ${same_hand} icons/same_same
overlay ${same_hemi} ${diff_hand} icons/same_diff
overlay ${diff_hemi} ${same_hand} icons/diff_same
overlay ${diff_hemi} ${diff_hand} icons/diff_diff
