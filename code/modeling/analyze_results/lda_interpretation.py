import pprint as pp
import pickle
import pandas as pd
import matplotlib as mpl
import numpy as np

from glob import glob
from pathlib import Path
from matplotlib import colors
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.inspection import DecisionBoundaryDisplay

# Plotting functions from here: https://scikit-learn.org/stable/auto_examples/classification/plot_lda_qda.html#sphx-glr-auto-examples-classification-plot-lda-qda-py

def plot_ellipse(mean, cov, color, ax):
    v, w = np.linalg.eigh(cov)
    u = w[0] / np.linalg.norm(w[0])
    angle = np.arctan(u[1] / u[0])
    angle = 180 * angle / np.pi  # convert to degrees
    # filled Gaussian at 2 standard deviation
    ell = mpl.patches.Ellipse(
        mean,
        2 * v[0] ** 0.5,
        2 * v[1] ** 0.5,
        angle=180 + angle,
        facecolor=color,
        edgecolor="black",
        linewidth=2,
    )
    ell.set_clip_box(ax.bbox)
    ell.set_alpha(0.4)
    ax.add_artist(ell)


def plot_result(estimator, X, y, ax):
    cmap = colors.ListedColormap(["tab:red", "tab:blue"])
    DecisionBoundaryDisplay.from_estimator(
        estimator,
        X,
        response_method="predict_proba",
        plot_method="pcolormesh",
        ax=ax,
        cmap="RdBu",
        alpha=0.3,
    )
    DecisionBoundaryDisplay.from_estimator(
        estimator,
        X,
        response_method="predict_proba",
        plot_method="contour",
        ax=ax,
        alpha=1.0,
        levels=[0.5],
    )
    y_pred = estimator.predict(X)
    X_right, y_right = X[y == y_pred], y[y == y_pred]
    X_wrong, y_wrong = X[y != y_pred], y[y != y_pred]
    ax.scatter(X_right[:, 0], X_right[:, 1], c=y_right, s=20, cmap=cmap, alpha=0.5)
    ax.scatter(
        X_wrong[:, 0],
        X_wrong[:, 1],
        c=y_wrong,
        s=30,
        cmap=cmap,
        alpha=0.9,
        marker="x",
    )
    ax.scatter(
        estimator.means_[:, 0],
        estimator.means_[:, 1],
        c="yellow",
        s=200,
        marker="*",
        edgecolor="black",
    )

    if isinstance(estimator, LinearDiscriminantAnalysis):
        covariance = [estimator.covariance_] * 2
    else:
        covariance = estimator.covariance_
    plot_ellipse(estimator.means_[0], covariance[0], "tab:red", ax)
    plot_ellipse(estimator.means_[1], covariance[1], "tab:blue", ax)

    ax.set_box_aspect(1)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set(xticks=[], yticks=[])

# My work

test_groups = ["A", "B", "C", "D", "E"]

# Set up paths

modeling = Path("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting/code/modeling/")
results_base = modeling / "results" / "base"

# Get HC to predict

with open(modeling / "hemiconnectome.pickle", 'rb') as f:
    hc = pickle.load(f)

# Get LDAs by looping over test group
ldas = []
for i, j in enumerate(test_groups):

    path = results_base / \
        f"method-lda_outcome-hemi_test-{j}_hands-all_data-base.pickle"

    with open(path, 'rb') as f:
        data = pickle.load(f)

    ldas += [data]

# Predict
print("Starting prediction ...")

labels_to_drop = ["sub", "handedness", "group", "gender", "age", "age_group",
                  "hemi", "class", "EHI"]

all_results = pd.DataFrame()
for i, j in enumerate(test_groups):

    test_data = hc[hc["group"] == j]

    identifiers = test_data[["sub", "group", "class"]].reset_index(drop=True)
    test_data = test_data.drop(labels_to_drop, axis=1, errors="ignore")

    results = pd.DataFrame(ldas[i].transform(test_data))
    results.columns = ["LD1", "LD2", "LD3"]
    results.reset_index(drop=True)

    results2 = pd.concat([identifiers, results], axis=1)

    all_results = pd.concat([all_results, results2])

all_results.reset_index(drop=True)

all_results.to_csv(modeling / "analyze_results" / "all_hemis_LDs.csv")
