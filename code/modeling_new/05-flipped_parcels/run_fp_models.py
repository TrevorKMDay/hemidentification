import pickle as pkl
import os
import sys
import pandas as pd

from pathlib import Path

sys.path.append("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting/code/" + 
                "modeling_new")

from run_model_kfold import run_my_kfold
from sklearn.feature_selection import SelectKBest, f_classif

home = Path("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting/code/",
            "modeling_new/05-flipped_parcels")

os.chdir(home)

def show_mcc(result):

    _, mcc, ci = result
    print(f"MCC: {mcc:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]")

# Load hemiconnectome

ll_file = home / ".." / "inputs"  / "hemiconnectome_LL.pickle"
with open(ll_file, "rb") as f:
    ll0 = pkl.load(f)

ll = ("hc.fpLL.RHD", ll0[ll0["handedness"] == "righty"])
ll_lefty = ("hc.fpLL.LHD", ll0[ll0["handedness"] == "lefty"])

rr_file = home / ".." / "inputs"  / "hemiconnectome_RR.pickle"
with open(rr_file, "rb") as f:
    rr0 = pkl.load(f)

rr = ("hc.fpRR.RHD", rr0[rr0["handedness"] == "righty"])
rr_lefty = ("hc.fpRR.LHD", rr0[rr0["handedness"] == "lefty"])

# Run models 

# Get columns that don't start with cx_
non_pred_cols = [x for x in ll0.columns if "cx_" not in x]

for i in ["lda"]: #, "nn", "svc"]:

    for j in [ll, rr, ll_lefty, rr_lefty]:
        
        os.makedirs(home / "results", exist_ok=True)
        
        run_my_kfold(j, 
                     outcome_cols="hemi", 
                     id_cols=non_pred_cols,
                     classifier=i,
                     output_dir=home / "results", 
                     bootstrap_mcc=True)

        y = j[1]["hemi"]
        x = j[1].drop(non_pred_cols, axis=1)

        selK = SelectKBest(f_classif, k=16).fit(x, y)
        
        selK_results = pd.DataFrame({"feature": selK.feature_names_in_,
                                     "score": selK.scores_})

        selK_results.to_csv(f"data-{j[0]}_selKbestscores.csv", index=False)