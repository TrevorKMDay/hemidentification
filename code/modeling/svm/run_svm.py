import pandas as pd
import pickle
import datetime as dt
from pathlib import Path

from sklearn import svm, metrics

# Setup =====

home = Path("/Users/tkmd/Google Drive/My Drive/Projects/"
            "hemisphere_fingerprinting")

# Load organized data from file
with open(home / "hemiconnectome.pickle", "rb") as f:
    data = pickle.load(f)

# print(data)

# Run model =====

# for test_group in ["A", "B", "C", "D", "E"]:

test_group = "A"

for k in ["linear", "poly", "rbf", "sigmoid"]:

    for w in [None, "balanced"]:

        print(f"{k} {w}")

        test_set = data[data["group"] == test_group]
        train_set = data[data["group"] != test_group]

        # Cut in half because there's two rows (LH/RH) per participant
        print(f"Sets of size test ({test_group}): {len(test_set) / 2 } and "
            f"training (!{test_group}): {len(train_set) / 2}\n")

        # lefty/righty L/R hemispheres
        train_labels = train_set["hemi"]
        # train_labels = pd.DataFrame({"hand": train_set["handedness"],
        #                              "hemi": train_set["hemi"]})
        training_data = train_set.drop(["sub", "handedness", "group", "gender",
                                        "age", "hemi", "class"],
                                    axis = 1)


        test_labels = test_set["hemi"]
        test_data = test_set.drop(["sub", "handedness", "group", "gender",
                                        "age", "hemi", "class"],
                                    axis = 1)

        # print(training_data.shape)
        # print(train_labels.shape)

        # class_weight='balanced' weights classes inversely proportional to their
        #   size in the data
        clf = svm.SVC(kernel=k, class_weight=w)

        print(f"Starting at: {dt.datetime.now()}")

        clf.fit(training_data, train_labels)
        print(clf.class_weight_)

        print(f"Finished at: {dt.datetime.now()}")

        test_result = clf.predict(test_data)

        # print(f"Accuracy: {clf.score(test_data, test_labels)}")

        print(metrics.confusion_matrix(test_labels, test_result))

        mcc = metrics.matthews_corrcoef(test_labels, test_result)
        print(f"MCC: {round(mcc, 3)}")

        all_results = pd.DataFrame({"sub": test_set["sub"],
                                    "gt": test_labels,
                                    "y_pred": test_result})

        wt = "None" if w is None else w
        all_results.to_csv(f"svcresults_kernel-{k}_wt-{wt}.csv")

        # break