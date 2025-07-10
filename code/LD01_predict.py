import sklearn
import glob as g
import pickle
import pandas as pd
import numpy as np

from pathlib import Path

modeling = Path("modeling")

# Get files
files = g.glob("modeling/results/base/method-lda_outcome-hand_*.pickle")
files.sort()

# Load models
models = []

for f in files:
    with open(f, 'rb') as p:
        models += [pickle.load(p)]

# Load data

with open(modeling / "inputs" / "hemiconnectome.pickle", "rb") as p:
    hc = pickle.load(p)

# Predict each fold

labels_to_drop = ["sub", "handedness", "handedness2", "group", "gender",
                  "age_group", "hemi", "class", "class2", "EHI"]

for i, group in enumerate(["A", "B", "C", "D", "E"]):

    # Get the fold from the data
    test = hc[hc["group"] == group]

    # Separate into two DFS
    labels = test[labels_to_drop]
    test = test.drop(labels_to_drop, axis = 1, errors = "ignore")

    # Do the prediction
    the_model = models[i]
    results = the_model.predict(test)
    xform = pd.DataFrame(the_model.transform(test))

    # Assign the result column to the labels for later processing
    labels["predicted"] = results

    # Combine labels with xform
    xform2 = pd.concat([labels.reset_index(drop=True), xform], axis=1)

    # feature weights
    names = pd.DataFrame(the_model.feature_names_in_, columns=["feature"])
    scalings = pd.DataFrame(the_model.scalings_,
                            columns=["scale_LD1", "scale_LD2", "scale_LD3"])
    coefs = pd.DataFrame(np.transpose(the_model.coef_),
                         columns=the_model.classes_,
                         index=the_model.feature_names_in_)


    features = pd.concat([names.reset_index(drop=True), scalings],
                         axis=1)

    print(features)

    # Write stuff out

    filename = f"method-lda_test-{group}"

    labels.to_csv(modeling / "results" / "base" / (filename + "_predicted.csv"),
                  index=False)

    xform2.to_csv(modeling / "results" / "base" / (filename + "_xform.csv"),
                  index=False)

    scalings.to_csv(modeling / "results" / "base" /
                     (filename + "_scalings.csv"),
                     index=False)

    coefs.to_csv(modeling / "results" / "base" / (filename + "_coefs.csv"))