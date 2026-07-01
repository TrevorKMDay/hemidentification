
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
            "modeling_new/06-handedness")

os.chdir(home)

def show_mcc(result):

    _, mcc, ci = result
    print(f"MCC: {mcc:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]")

# Load hemiconnectome

hc_file = home / ".." / "inputs"  / "hc_hpZ.pickle"
with open(hc_file, "rb") as f:
    hc0 = pkl.load(f)

hc = ("hc", hc0)

# Run models 

# Get columns that don't start with cx_
non_pred_cols = [x for x in hc0.columns if "cx_" not in x]

# Run the all-in model and keep the scalings 
run_my_kfold(hc, 
             outcome_cols="class", 
             id_cols=non_pred_cols,
             group_col=None,
             classifier="lda",
             output_dir=home / "all_in", 
             bootstrap_mcc=True)
