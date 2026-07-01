import pickle as pkl

from run_model_kfold import create_null_distribution
from pathlib import Path

home = Path("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting/code/",
            "modeling_new/07-bootstrapping")

hc_file = home / ".." / "inputs"  / "hc_hpZ.pickle"
with open(hc_file, "rb") as f:
    hc0 = pkl.load(f)

hc = ("hc", hc0)

# Get columns that don't start with cx_
non_pred_cols = [x for x in hc0.columns if "cx_" not in x]

create_null_distribution(hc, 
                         id_cols=non_pred_cols,
                         cols_to_shuffle=["hemi"], 
                         output_dir=home / "bs_results",
                         n=100)

