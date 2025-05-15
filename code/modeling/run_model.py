import pandas as pd
import pickle
import datetime as dt
import pprint as pp
import os
import re

from pathlib import Path

import sklearn
from sklearn import metrics

import argparse as ap

import sklearn.linear_model
import sklearn.neural_network

from random import sample

parser = ap.ArgumentParser()

parser.add_argument("predict",
                    help="Column to predict")

parser.add_argument("data",
                    help="Pickle file with data.")

parser.add_argument("test_group", choices=["A", "B", "C", "D", "E", "within"],
                    help="Which test group to hold out.")

parser.add_argument("method", choices=["sv", "nn", "lda"],
                    help="Which test group to hold out.")

parser.add_argument("--hands", choices=["lefty", "righty"],
                    help="Select which handedness group to KEEP.")

parser.add_argument("--hemis", choices=["LH", "RH"],
                    help="Select which hemispheres to KEEP")

parser.add_argument("--output_dir", "-o",
                    help="Dir to save file to. Default: .")

parser.add_argument("--output_name",
                    help="Full path to save to.")

parser.add_argument("--bootstrap", "-b", type=int, metavar="N",
                    help="From the data provided, sample N columns.")

parser.add_argument("--force", action="store_true")

args = parser.parse_args()

outcome = args.predict
input_data = args.data
test_group = args.test_group
method = args.method
hands = args.hands
hemis = args.hemis
bootstrap = args.bootstrap

opath = Path(args.output_dir) if args.output_dir is not None else None
oname = args.output_name

# Setup =====

home = Path("/Users/tkmd/Google Drive/My Drive/Projects/"
            "hemisphere_fingerprinting")

hand_suffix = f"_hands-{hands}" if hands is not None else ""
hemi_suffix = f"_hemi-{hemis}" if hemis is not None else ""
bs_suffix = f"_bs-{bootstrap}" if bootstrap is not None else ""

# Check output ====

if oname is None:

    # If no output path, use cwd
    opath = Path(".") if opath is None else opath

    output_name = opath / \
                (f"method-{method}_outcome-{outcome}_" + \
                f"group-{test_group}" + hand_suffix + hemi_suffix + \
                bs_suffix + ".csv")

else:

    output_name = oname

if os.path.exists(output_name) and not args.force:
    print(f"Output file {output_name} exists!")
    raise SystemExit

# Load data ====

# Load organized data from file
with open(input_data, "rb") as f:
    data = pickle.load(f)

print(f"Data loaded from {input_data}")

# Check outcome =====

if outcome in data.columns:

    outcome_type = data.dtypes[outcome]
    print(f"Requested outcome data type is: {outcome_type}")
else:

    print(f"Requested outcome '{outcome}' is missing from input data")
    raise SystemExit

if outcome_type == "object":

    if method == "sv":
        clf = sklearn.svm.SVC(kernel="linear", class_weight="balanced")
    elif method == "nn":
        clf = sklearn.neural_network.MLPClassifier(random_state=1, max_iter=300)
    elif method == "lda":
        clf = sklearn.discriminant_analysis.LinearDiscriminantAnalysis()

elif outcome_type == "float64":

    if method == "sv":
        clf = sklearn.svm.SVR(kernel="linear")
    elif method == "nn":
        clf = sklearn.neural_network.MLPRegressor(random_state=1, max_iter=300)
    elif method == "lda":
        print("Method lda not available for continuous outcome")
        raise SystemExit

print()
print("Starting modeling run ...")
print(f"Using method {method}!")


# Run model =====

if hands is not None:
    if "handedness" in data.columns:
        data = data[data["handedness"] == hands]
    else:
        print("Filter by handedness was asked for, but no 'handedness' column")
        raise SystemExit

if hemis is not None:
    if "hemi" in data.columns:
        data = data[data["hemi"] == hemis]
    else:
        print("Filter by hemisphere was asked for, but no 'hemi' column")
        raise SystemExit

# Do the modeling

if test_group != "within":

    test1 = data[data["group"] == test_group]
    train1 = data[data["group"] != test_group]

else:

    # If we are doing all-to-all, just use both as train and test
    test1 = train1 = data

print()
print("Train/test info:")

print("  Training data:")

if outcome_type == "object":
    print(train1[outcome].value_counts())

labels_to_drop = ["sub", "handedness", "handedness2", "group", "gender", "age",
                  "age_group", "hemi", "class", "class2", "EHI"]

# Create groups
train_labels = train1[outcome]
train_data = train1.drop(labels_to_drop, axis = 1, errors = "ignore")

exit

# Drop duplicates from oversampled data
test1.drop_duplicates(inplace=True)

test_labels = test1[outcome]
test_data = test1.drop(labels_to_drop, axis = 1, errors = "ignore")

label_check1 = [x for x in train_data.columns if not re.fullmatch(".*_.*", x)]
label_check2 = [x for x in test_data.columns if not re.fullmatch(".*_.*", x)]

print("\n  Possible bad columns:")
pp.pprint([label_check1, label_check2], compact=True)

print("  Test data:")

if outcome_type == "object":
    print(test_labels.value_counts())

if bootstrap is not None:

    cols_to_keep = sample(train_data.columns.tolist(), bootstrap)

    train_data = train_data[cols_to_keep]
    test_data = test_data[cols_to_keep]

    print(f"    Reduced train/test data to {bootstrap} columns.")

print(f"  Sets of size train: ({test_group}): "
      f"{len(train_data)}x{len(train_data.columns)} and test "
      f"(!{test_group}): {len(test_data)}x{len(test_data.columns)}\n")

#
# Start running model
#

print()
print(f"Starting at: {dt.datetime.now()}")

clf.fit(train_data, train_labels)

print(f"Finished at: {dt.datetime.now()}")
print()

# Write out model

output_pickle = str(output_name).replace("csv", "pickle")
with open(output_pickle, 'wb') as f:
    pickle.dump(clf, f)

# Summary values

print("Doing prediction ...")
test_result = clf.predict(test_data)

if outcome_type == "object":

    # To do: Equivalent for float

    acc = clf.score(test_data, test_labels)
    print(f"Accuracy: {round(acc, 3)}")

    c = metrics.confusion_matrix(test_labels, test_result)
    print(c)

    mcc = metrics.matthews_corrcoef(test_labels, test_result)
    print(f"MCC: {round(mcc, 3)}")

# Summarize results

results = pd.DataFrame({"sub": test1["sub"],
                        "ground": test_labels,
                        "predicted": test_result})

results.to_csv(output_name)

print()
print(f"Saved results to {output_name}")
print("Done!")