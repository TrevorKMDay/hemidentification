# THE PURPOSE OF THIS SCRIPT IS TO RUN THE SMALL RIGHT-HANDED MODELS TO COMPARE
# WITH THE SMALL LEFT-HANDED MODELS.
# THE CREATION OF THE RIGHT-HANDED GROUPS IS DONE IN R.

import pickle as pkl
import os
import sys
import pandas as pd

from pathlib import Path

sys.path.append("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting/code/" + 
                "modeling_new")
from run_model_kfold import run_my_kfold

home = Path("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting/code/",
            "modeling_new/03-left_handed")

os.chdir(home)

# Load sample determinations

subs = pd.read_csv("../righty_subsample_subs.csv")

random_subs = [s for s, x in zip(subs["sub"], subs["random"]) if x]
paired_subs = [s for s, x in zip(subs["sub"], subs["paired"]) if x]

# Normalized within hemisphere

hchpZ_file = home / ".." / "inputs" / "hc_hpZ.pickle"
with open(hchpZ_file, "rb") as f:
    hchpZ0 = pkl.load(f)

hchpZ0_random = hchpZ0[hchpZ0["sub"].isin(random_subs)]
hchpZ_random = ("hc.hpZ.RHDrdm", hchpZ0_random)

hchpZ0_paired = hchpZ0[hchpZ0["sub"].isin(paired_subs)]
hchpZ_paired = ("hc.hpZ.RHDprd", hchpZ0_paired)

lang_file = home / ".." / "inputs" / "hemiconnectome_language_hpZ.pickle"
with open(lang_file, "rb") as f:
    lang0 = pkl.load(f)

langZ0_random = lang0[lang0["sub"].isin(random_subs)]
langZ_random = ("lang.hpZ.RHDrdm", langZ0_random)

langZ0_paired = lang0[lang0["sub"].isin(paired_subs)]
langZ_paired = ("lang.hpZ.RHDprd", langZ0_paired)

# One hemi per

hchpZ_1p_random = ("hc.hpZ.1p.RHDrdm", 
                    hchpZ0_random[hchpZ0_random["hemi"] == hchpZ0_random["hemi_fold"]])

hchpZ_1p_paired = ("hc.hpZ.1p.RHDprd", 
                    hchpZ0_paired[hchpZ0_paired["hemi"] == hchpZ0_paired["hemi_fold"]])

langZ_1p_random = ("lang.hpZ.1p.RHDrdm", 
                    langZ0_random[langZ0_random["hemi"] == langZ0_random["hemi_fold"]])

langZ_1p_paired = ("lang.hpZ.1p.RHDprd", 
                    langZ0_paired[langZ0_paired["hemi"] == langZ0_paired["hemi_fold"]])

# Run models 

# Get columns that don't start with cx_
non_pred_cols = [x for x in hchpZ0.columns if "cx_" not in x]

for i in ["lda", "nn", "svc"]:

    for j in [hchpZ_random, hchpZ_paired, hchpZ_1p_random, hchpZ_1p_paired,
                langZ_random, langZ_paired, langZ_1p_random, langZ_1p_paired]:
        
        os.makedirs(home / "results" / "righty_pair", exist_ok=True)
        run_my_kfold(j, 
                    outcome_cols="hemi", 
                    id_cols=non_pred_cols,
                    classifier=i,
                    output_dir=home / "results" / "righty_pair", 
                    bootstrap_mcc=True)

