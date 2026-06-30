import pickle as pkl

from run_model_kfold import run_my_kfold, bootstrap_mcc_value
from pathlib import Path

home = Path("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting/code/",
            "modeling_new")

def show_mcc(result):

    _, mcc, ci = result
    print(f"MCC: {mcc:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]")


non_pred_cols = ("sub", "hemi_fold", "age_group", "hemi", "gender",
                 "handedness", "EHI", "group")

# Completely random data 

with open(home / "inputs" / "completely_random.pickle", 'rb') as f:
    crandom_data = pkl.load(f)

crandom = ("random", crandom_data)
crandom_result = run_my_kfold(crandom, outcome_cols="hemi", 
                                id_cols=non_pred_cols,
                                output_dir=home / "results", 
                                bootstrap_mcc=True)

# Shuffled within hemispheres

with open(home / "inputs" / "shuffled_labels" / "shuffled_separately.pickle", 'rb') as f:
    hemi_random = ("hemi_random", pkl.load(f))

crandom_result = run_my_kfold(hemi_random, outcome_cols="hemi", 
                                id_cols=non_pred_cols,
                                output_dir=home / "results", 
                                bootstrap_mcc=True)

# Scalars 

with open(home / "inputs" / "summary_values.pickle", 'rb') as f:
    summary_values = ("summary", pkl.load(f))

run_my_kfold(summary_values,
             outcome_cols="hemi", id_cols=["sub"],
             output_dir=home / "results" / "shuffled_labels", 
             bootstrap_mcc=True)

vars = summary_values[1].columns[3:7]

for v in vars:

    print()
    print(f"{v}:")
    temp = (f"summary_{v}", summary_values[1][["group", "sub", "hemi", v]])

    run_my_kfold(temp,
                 outcome_cols="hemi", id_cols=["sub"],
                 output_dir=home / "results" / "shuffled_labels", 
                 bootstrap_mcc=True)

# Load data into a list of data-label pairs

data = []

for i in range(1, 11):

    bname = f"shuffle{i:02}"
    file = home / "inputs" / "shuffled_labels" / f"{bname}.pickle"
    
    with open(file, 'rb') as f:
        temp = pkl.load(f)
    
    data_label = [(bname, temp)]
    data += data_label

# Set up results

results = []

for i in range(0, 11):

    r = run_my_kfold(data[i], outcome_cols="hemi", 
                        id_cols=non_pred_cols,
                        output_dir=home / "results" / "shuffled_labels", 
                        bootstrap_mcc=True)

    results += [r]
