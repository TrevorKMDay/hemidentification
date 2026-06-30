
# DIMENSIONALITIY REDUCTION FOR HCP DATA

# import sklearn
import pickle as pkl
import pandas as pd

from pathlib import Path

# Tools from scikit learn

# from sklearn import metrics
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.feature_selection import  SelectKBest, f_classif

# Set up directories ====
home = Path("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting")
data_file = home / "code" / "modeling" / "inputs" / "hemiconnectome.pickle"
output_dir = home / "code" / "modeling" / "dim_reduce"

# Load the data ====

print(f"Loading data from {data_file}")
with open(data_file, "rb") as f:
    data = pd.DataFrame(pkl.load(f))

# How to reduce columns

def sel_kbest(k, train, test, y_name, verbose=False):

    print(f"k: {k}")

    train_x, train_y = train
    test_x, test_y = test

    selector = SelectKBest(f_classif, k=k)
    train_new = pd.DataFrame(selector.fit_transform(train_x, y=train_y))

    # Get the names of the features that were kept by the selector
    feature_mask = selector.get_support()
    selected_names = train_x.columns[feature_mask]
    train_new.columns = selected_names

    names_ordered = [val for _, val in
                    sorted(zip(selector.scores_, train_x.columns))]

    scores_ordered = sorted(selector.scores_)

    if verbose:
        for name, score in zip(names_ordered, scores_ordered):
            print(f"Feature: {name:<20} Score: {score:.2f}")

        print(f"\nChosen Features: {list(selected_names)}")

    # Run the LDA and report the within-sample data
    clf = LinearDiscriminantAnalysis()
    model = clf.fit(train_new, train_y)
    acc = model.score(train_new, train_y)

    print(f"Self accuracy ({y_name}): {round(acc, 3)}")

    test_new = test_x[selected_names]

    # Do the out-of-group prediction
    test_y["predicted"] = model.predict(test_new)
    acc2 = (test_y[y_name] == test_y['predicted']).mean()
    print(f"Test accuracy ({y_name}) {round(acc2, 3)}")

    return([model, test_y])

# Five-fold

test_groups = ["A", "B", "C", "D", "E"]
models_hemi = models_class = list()
pred_hemi = pred_class = pd.DataFrame()

for test_group in test_groups:

    test1 = data[data["group"] == test_group]
    train1 = data[data["group"] != test_group]

    id = train1[["sub", "hemi", "handedness"]]
    hemi_y = train1["hemi"].tolist()
    class_y = train1["class"].tolist()

    print(f"Testing on {test_group}, training on the rest ...")

    # Training data
    features = train1.drop(["sub", "gender", "age", "age_group", "hemi",
                            "handedness", "handedness2", "class", "class2",
                            "group", "EHI"],
                           axis=1, errors="ignore")

    # Subset the feature columns
    test_x = test1[train1.columns]
    test_y = test1[["sub", "hemi"]]
    test_y2 = test1[["sub", "class"]]

    for k in [2**x for x in range(10)]:

        # Test how few features work on everyone

        m1, y1 = sel_kbest(k=k, train=[features, hemi_y],
                           test=[test_x, test_y],
                           y_name = "hemi")

        # Collect up results
        models_hemi.append([[test_group, k, m1]])

        y1["group"] = test_group
        y1["k"] = k
        pred_hemi = pd.concat([pred_hemi, y1])

        m2, y2 = sel_kbest(k=k, train=[features, class_y],
                           test=[test_x, test_y2],
                           y_name = "class")

        models_class.append([[test_group, k, m2]])

        y2["group"] = test_group
        y2["k"] = k
        pred_class = pd.concat([pred_class, y2])

pred_hemi.to_csv(output_dir / "data-all_outcome-hemi_results.csv")
pred_class.to_csv(output_dir / "data-all_outcome-class_results.csv")

modelsL = list()
pred_hemiL = pd.DataFrame()


for test_group in test_groups:

    # Do this just for lefties

    testL = data[(data["group"] == test_group) & (data["handedness"] == "lefty")]
    trainL = data[(data["group"] != test_group) & (data["handedness"] == "lefty")]

    id = trainL[["sub", "hemi", "handedness"]]
    hemi_y = trainL["hemi"].tolist()

    print(f"Testing on {test_group}, training on the rest ...")

    # Training data
    features = trainL.drop(["sub", "gender", "age", "age_group", "hemi",
                            "handedness", "handedness2", "class", "class2",
                            "group", "EHI"],
                           axis=1, errors="ignore")

    # Subset the feature columns
    testLx = testL[trainL.columns]
    testLy = testL[["sub", "hemi"]]

    for k in [2**x for x in range(13)]:

        # Test how few features work on everyone

        m1, y1 = sel_kbest(k=k, train=[features, hemi_y],
                           test=[testLx, testLy],
                           y_name = "hemi")

        # Collect up results
        modelsL.append([[test_group, k, m1]])

        y1["group"] = test_group
        y1["k"] = k
        pred_hemiL = pd.concat([pred_hemiL, y1])
