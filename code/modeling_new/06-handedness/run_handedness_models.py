
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

hco_file = home / ".." / "inputs"  / "hc_oversampled.pickle"
with open(hco_file, "rb") as f:
    hco0 = pkl.load(f)

hco = ("hco", hco0)

trans_file = home / ".." / "inputs"  / "transconnectome.pickle"
with open(trans_file, "rb") as f:
    trans0 = pkl.load(f)

trans = ("trans", trans0)

fc_file = home / ".." / "inputs"  / "connectome.pickle"
with open(fc_file, "rb") as f:
    fc0 = pkl.load(f)

fc = ("fc", fc0)

lg_file = home / ".." / "inputs"  / "hemiconnectome_language_hpZ.pickle"
with open(lg_file, "rb") as f:
    lang0 = pkl.load(f)

lang = ("lang", lang0)

# Run models 

# Get columns that don't start with cx_
non_pred_cols = [x for x in hc0.columns if "cx_" not in x]

for i in ["lda", "nn", "svc"]:

    for j in [hc, hco, lang]:
        
        os.makedirs(home / "results", exist_ok=True)
        
        run_my_kfold(j, 
                     outcome_cols="class", 
                     id_cols=non_pred_cols,
                     classifier=i,
                     output_dir=home / "results", 
                     bootstrap_mcc=True)

    run_my_kfold(hc, 
                 outcome_cols=["class", "gender"], 
                 id_cols=non_pred_cols,
                 classifier=i,
                 output_dir=home / "results", 
                 bootstrap_mcc=True)

    for j in [fc, trans]:
        
        os.makedirs(home / "results", exist_ok=True)
        
        run_my_kfold(j, 
                     outcome_cols="handedness", 
                     id_cols=non_pred_cols,
                     classifier=i,
                     output_dir=home / "results", 
                     bootstrap_mcc=True)
                     
        run_my_kfold(j, 
                     outcome_cols=["handedness", "gender"],
                     id_cols=non_pred_cols,
                     classifier=i,
                     output_dir=home / "results", 
                     bootstrap_mcc=True)

