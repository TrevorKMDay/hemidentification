import pickle as pkl

from run_model_kfold import run_my_kfold, bootstrap_mcc_value
from pathlib import Path

home = Path("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting/code/",
            "modeling_new")

def show_mcc(result):

    _, mcc, ci = result
    print(f"MCC: {mcc:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]")

# Load hemiconnectome

hc_file = home / "inputs" / "hemiconnectome.pickle"
with open(hc_file, "rb") as f:
    hc0 = pkl.load(f)

hemiconnectome = ("hemiconnectome", hc0)
hemi_lhonly = ("hc_lhonly", hc0[hc0["hemi"] == "LH"])
hemi_rhonly = ("hc_rhonly", hc0[hc0["hemi"] == "RH"])

non_pred_cols = ("sub", "hemi_fold", "age_group", "hemi", "gender",
                 "handedness", "class", "EHI", "group")

# Basic gender models ====

gender_bi = run_my_kfold(hemiconnectome, outcome_cols="gender", 
                         id_cols=non_pred_cols,
                         output_dir=home / "results/", bootstrap_mcc=True)

gender_lh = run_my_kfold(hemi_lhonly, outcome_cols="gender", 
                         id_cols=non_pred_cols,
                         output_dir=home / "results/", bootstrap_mcc=True)

gender_rh = run_my_kfold(hemi_rhonly, outcome_cols="gender", 
                         id_cols=non_pred_cols,
                         output_dir=home / "results/", bootstrap_mcc=True)

# Basic age models ====

age_bi = run_my_kfold(hemiconnectome, outcome_cols="age_group", 
                         id_cols=non_pred_cols,
                         output_dir=home / "results/", bootstrap_mcc=True)

age_lh = run_my_kfold(hemi_lhonly, outcome_cols="age_group", 
                         id_cols=non_pred_cols,
                         output_dir=home / "results/", bootstrap_mcc=True)

age_rh = run_my_kfold(hemi_rhonly, outcome_cols="age_group", 
                         id_cols=non_pred_cols,
                         output_dir=home / "results/", bootstrap_mcc=True)

# Basic hand models ====

hand_bi = run_my_kfold(hemiconnectome, outcome_cols="handedness", 
                         id_cols=non_pred_cols,
                         output_dir=home / "results/", bootstrap_mcc=True)

hand_lh = run_my_kfold(hemi_lhonly, outcome_cols="handedness", 
                         id_cols=non_pred_cols,
                         output_dir=home / "results/", bootstrap_mcc=True)

hand_rh = run_my_kfold(hemi_rhonly, outcome_cols="handedness", 
                         id_cols=non_pred_cols,
                         output_dir=home / "results/", bootstrap_mcc=True)

# Basic hemisphere model ====

hemi = run_my_kfold(hemiconnectome, outcome_cols="hemi", 
                         id_cols=non_pred_cols,
                         output_dir=home / "results/", bootstrap_mcc=True)

# Crossing demographics -----

gender_age = run_my_kfold(hemiconnectome, outcome_cols=["gender", "age_group"],
                         id_cols=non_pred_cols,
                         output_dir=home / "results/", bootstrap_mcc=True)

gender_hand = run_my_kfold(hemiconnectome, 
                           outcome_cols=["gender", "handedness"],
                           id_cols=non_pred_cols,
                           output_dir=home / "results/", bootstrap_mcc=True)

age_hand = run_my_kfold(hemiconnectome, 
                        outcome_cols=["age_group", "handedness"],
                        id_cols=non_pred_cols,
                        output_dir=home / "results/", bootstrap_mcc=True)

all = run_my_kfold(hemiconnectome, 
                        outcome_cols=["gender", "age_group", "handedness"],
                        id_cols=non_pred_cols,
                        output_dir=home / "results/", bootstrap_mcc=True)

# Crossing hemisphere =====

gender_hemi = run_my_kfold(hemiconnectome, outcome_cols=["gender", "hemi"],
                           id_cols=non_pred_cols,
                           output_dir=home / "results/", bootstrap_mcc=True)

age_hemi = run_my_kfold(hemiconnectome, 
                           outcome_cols=["age_group", "hemi"],
                           id_cols=non_pred_cols,
                           output_dir=home / "results/", bootstrap_mcc=True)

hand_hemi = run_my_kfold(hemiconnectome, 
                        outcome_cols=["handedness", "hemi"],
                        id_cols=non_pred_cols,
                        output_dir=home / "results/", bootstrap_mcc=True)

# Check for leakage ----

one_hemi_each = ("one_hemi_each",
    hemiconnectome[1][hemiconnectome[1]["hemi"] == hemiconnectome[1]["hemi_fold"]])

hemi_oneeach =run_my_kfold(one_hemi_each, 
                            outcome_cols=["hemi"],
                            id_cols=non_pred_cols,
                            output_dir=home / "results/", bootstrap_mcc=True)
