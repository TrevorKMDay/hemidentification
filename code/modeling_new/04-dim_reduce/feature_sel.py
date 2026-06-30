# from os import F_OK
import pickle as pkl 
import pandas as pd
import math
import matplotlib.pyplot as plt

# This file is symlinked to this directory, probably hacky but whatever
from run_model_kfold import run_my_kfold, bootstrap_mcc_value
from pathlib import Path
from sklearn.feature_selection import SelectKBest, f_classif


home = Path("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting/code/"
            "modeling_new/04-dim_reduce")

# Matched samples 
subs = pd.read_csv(home / ".." / "righty_subsample_subs.csv")
paired_subs = [s for s, x in zip(subs["sub"], subs["paired"]) if x]

hc_file = home / ".." / "inputs" / "hemiconnectome.pickle"
with open(hc_file, "rb") as f:
    hc0 = pkl.load(f)

hemiconnectome = ("hemiconnectome", hc0)
hc_righties = ("hcrhd", hc0[hc0["handedness"] == "righty"])
hc_rhdpaired = ("hcrhd.prd", hc0[hc0["sub"].isin(paired_subs)])
hc_lefties = ("hclhd", hc0[hc0["handedness"] == "lefty"])

# Common rule of thumb
rule_10x = [int(math.floor(x[1].shape[0] / 10)) 
            for x in [hc_righties, hc_lefties]]

# For highly-correlated features
rule_sqrt = [math.ceil(math.sqrt(x[1].shape[0]))
             for x in [hc_righties, hc_lefties]]

# Does hemidentification work with a lot fewer features?
non_pred_cols = ["sub", "hemi_fold", "age_group", "hemi", "gender",
                 "handedness", "class", "EHI", "group"]

hc_id = hemiconnectome[1][non_pred_cols]
hc_x = hemiconnectome[1].drop(non_pred_cols, axis=1)
hc_y = hc_id["hemi"]

# Feature selection over sizes

# Do 2^x but also the additional sizes
sizes = [2 ** x for x in range(14)] + rule_10x + rule_sqrt
sizes.sort()

result_list = []
for dataset in [hc_lefties, hc_righties, hc_rhdpaired]:

    id = dataset[1][non_pred_cols]
    x = dataset[1].drop(non_pred_cols, axis=1)
    y = id["hemi"]

    for size in sizes:

        print(f"{dataset[0]} {size}")

        selK = SelectKBest(f_classif, k=size).fit(x, y)
        temp_data = (f"{dataset[0]}_k-{size}", 
                    pd.concat([id, x[selK.get_feature_names_out()]],
                              axis=1))

        results = run_my_kfold(temp_data, outcome_cols="hemi",
                                id_cols=non_pred_cols, classifier="lda",
                                output_dir=home / "outputs",
                                bootstrap_mcc=True, collapse_mcc=False)

        result_list += [[dataset[0], size, results]]

with open(home / "feature_sel_results.pickle", "wb") as f:
    pkl.dump(result_list, f)

# Get just the feature scores

top_lefties = SelectKBest(f_classif, k=16).fit(hc_lefties[1].drop(non_pred_cols, axis=1),
                                                hc_lefties[1]["hemi"])

top_lefties_table = pd.DataFrame({"feature": top_lefties.feature_names_in_,
                                  "lefty_score": top_lefties.scores_})

top_lefties_table.to_csv(home / "lefty_scores.csv", index=False)

top_righties = SelectKBest(f_classif, k=16).fit(hc_righties[1].drop(non_pred_cols, axis=1),
                                                hc_righties[1]["hemi"])

top_righties_table = pd.DataFrame({"feature": top_righties.feature_names_in_,
                                   "righty_score": top_righties.scores_})

top_righties_table.to_csv(home / "righty_scores.csv", index=False)

# Run the within-model prediction
results_list = []
for dataset in [hemiconnectome, hc_lefties, hc_righties, hc_rhdpaired]:

    id = dataset[1][non_pred_cols]
    x = dataset[1].drop(non_pred_cols, axis=1)
    y = id["hemi"]

    for size in sizes + [16110]:

        print(f"{dataset[0]} {size}")

        selK = SelectKBest(f_classif, k=size).fit(x, y)
        temp_data = (f"{dataset[0]}_k-{size}", 
                    pd.concat([id, x[selK.get_feature_names_out()]],
                              axis=1))

        results = run_my_kfold(temp_data, outcome_cols="hemi",
                                id_cols=non_pred_cols, classifier="lda",
                                output_dir=home / "all_in",
                                group_col=None,
                                bootstrap_mcc=True, collapse_mcc=True,
                                dump_scalings=True)

        result_list += [[dataset[0], size, results]]
