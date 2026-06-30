import pickle as pkl
import os
import sys

from pathlib import Path

sys.path.append("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting/code/" + 
                "modeling_new")
from run_model_kfold import run_my_kfold

home = Path("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting/code/",
            "modeling_new/03-right_handed")

os.chdir(home)

def show_mcc(result):

    _, mcc, ci = result
    print(f"MCC: {mcc:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]")

# Load hemiconnectome

hc_file = home / ".." / "inputs"  / "hemiconnectome.pickle"
with open(hc_file, "rb") as f:
    hc0 = pkl.load(f)

hc0 = hc0[hc0["handedness"] == "righty"]
hemiconnectome = ("hc.RHD", hc0)

# Normalized within hemisphere

hchpZ_file = home / ".." / "inputs" / "hc_hpZ.pickle"
with open(hchpZ_file, "rb") as f:
    hchpZ0 = pkl.load(f)

hchpZ0 = hchpZ0[hchpZ0["handedness"] == "righty"]
hchpZ = ("hc.hpZ.RHD", hchpZ0)

lang_file = home / ".." / "inputs" / "hemiconnectome_language_hpZ.pickle"
with open(lang_file, "rb") as f:
    lang0 = pkl.load(f)

lang0 = lang0[lang0["handedness"] == "righty"]
langhpZ = ("lang.hpZ.RHD", lang0)

# One hemi per

hc_1p = ("hc.1p.RHD", hc0[hc0["hemi"] == hc0["hemi_fold"]])
hchpZ_1p = ("hc.hpZ.1p.RHD", hchpZ0[hchpZ0["hemi"] == hchpZ0["hemi_fold"]])
langhpZ_1p = ("lang.hpZ.1p.RHD", lang0[lang0["hemi"] == lang0["hemi_fold"]])

# Run models 

# Get columns that don't start with cx_
non_pred_cols = [x for x in hc0.columns if "cx_" not in x]

for i in ["lda", "nn", "svc"]:

    for j in [hemiconnectome, hchpZ, hc_1p, hchpZ_1p, langhpZ, langhpZ_1p]:
        
        os.makedirs(home / "results" / "hemiconnectome", exist_ok=True)
        
        run_my_kfold(j, 
                    outcome_cols="hemi", 
                    id_cols=non_pred_cols,
                    classifier=i,
                    output_dir=home / "results" / "hemiconnectome", 
                    bootstrap_mcc=True)

