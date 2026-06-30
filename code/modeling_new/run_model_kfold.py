
from matplotlib.pyplot import axis
from IPython.core import error
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import copy
import sklearn
import os

from sklearn.metrics import matthews_corrcoef
from pathlib import Path

def bootstrap_mcc_value(y_true, y_pred, n_reps=500):

    df_orig = pd.DataFrame({"y_true": y_true, "y_pred": y_pred})

    # Sample n reps
    mccs = list()
    for i in range(n_reps):
        df_new = df_orig.sample(frac = 1, replace=True)
        mcc_new = matthews_corrcoef(y_true=df_new["y_true"], 
                                    y_pred=df_new["y_pred"])
        mccs += [mcc_new]

    mean_mcc = np.mean(mccs)
    mcc_quantiles = np.quantile(mccs, [0.025, 0.975])

    # Plot
    # fig, ax = plt.subplots()

    # ax.hist(mccs, bins=8, linewidth=0.5, edgecolor="white")
    # ax.set(xlim=(0, 1))
    # plt.show()

    return((mean_mcc, mcc_quantiles))

def run_my_kfold(data_tuple, outcome_cols, id_cols, classifier,
                 output_dir, group_col="group", bootstrap_mcc=False,
                 collapse_mcc=True, dump_scalings=False):

    """
    Run k-fold prediction on a data-frame.

    Args:

        data_tuple (tuple): A 2-tuple that contains the label and a data frame
            with at least the following columns.

        outcome_cols (string): A string or list of strings containing the 
            columns to predict.

        id_cols (string): A list of strings containing which columns are 
            identifiers and not predictors. 

        classifier (string): lda, nn, or svc

        output_dir (Path object): Where to write out the files. 

        group_col (string): Default: "group", the column to base the folds on.
            If None, does all-to-all prediction.

        bootstrap_mcc (logical): If False, the raw MCC is returned without 
            a bootstrapped 95% CI.

        collapse_mcc (logical): If True, only overall MCC is returned, and not 
            on a per-group_col basis. 

        dump_scalings (logical): If True, the scalings (per feature are saved
            to a CSV). 

    Returns:

        A 3-list containing a data frame of predicted outcomes, per-group MCCs
        and the overall MCC, and 95% CIs, if requested.


    """

    # Split data name and the actual data from the tuple
    dt = copy.deepcopy(data_tuple)
    data_name, data = dt

    print(f"\nStarting work on: {data_name}")

    # Replace a string with one a one-item list
    if isinstance(outcome_cols, str):
        outcome_cols = [outcome_cols]

    print(f"    Dump scalings: {dump_scalings}")

    # Remove the control columns from the ID columns
    id_cols_local = list(id_cols)
    ctrl_cols = outcome_cols + [group_col]
    
    removed_cols = list(set(ctrl_cols) - set(id_cols_local))
    id_cols_local = [x for x in id_cols_local if x not in ctrl_cols]

    if len(removed_cols) > 0:
        print(f"    Removed control columns {removed_cols} from ID columns")

    removed_cols2 = list(set(id_cols_local) - set(data.columns))
    id_cols_local = [x for x in id_cols_local if x in data.columns]

    if len(removed_cols2) > 0:
        print(f"    Removed columns {removed_cols2} from ID columns")

    print(f"    Final ID columns: {id_cols_local}")

    # Create new outcome column by joining names
    all_outcomes = '.'.join(outcome_cols)
    # Create an outcome column for crossed vars if necessary
    if len(outcome_cols) > 1:
        outcome_df = data[outcome_cols]
        data[all_outcomes] = outcome_df.apply('.'.join, axis=1)

    print(f"    Predicting groups: {set(data[all_outcomes])}")

    # Create empty dataframe to store results. We're just gonna inefficently
    # append. 
    results_df = pd.DataFrame(columns=id_cols_local + ["predicted"])
    group_mccs = pd.DataFrame(columns=["group", "mcc", "mcc_ci_lo", 
                                       "mcc_ci_hi"])

    # Save results
    output_name = Path(output_dir) / \
        f"data-{data_name}_clf-{classifier}_outcome-{all_outcomes}_results.csv"

    scalings_name = Path(output_dir) / \
        f"data-{data_name}_clf-{classifier}_outcome-{all_outcomes}_scalings.csv"

    if os.path.exists(output_name):
        print(f"Requested output file {output_name} already exists!")
        return()

    # Get the groups from the group column
    if group_col is not None:
        groups = list(set(data[group_col]))
        groups.sort()
    else:
        groups = [None]

    print("    Running hold outs ...")
    for hold_out in groups:

        print(f"        {hold_out} ", end="")

        if hold_out is not None:
            # Filter data down to non-held-out group or only the held-out group
            # and then drop the grouping column from both train and test
            train = data[data[group_col] != hold_out].drop([group_col], axis=1, 
                                                           errors="ignore")
            test = data[data[group_col] == hold_out].drop([group_col], axis=1, 
                                                          errors="ignore")
        else:
            # If doing 'None', it's test/train within
            train = test = data

        # Extract the ground truths
        train_y = train[all_outcomes].tolist()
        test_y = test[all_outcomes].reset_index(drop=True)

        # Extract the columns to re-add to the output
        # train_id = train[id_cols_local]
        # CURRENT PROBLEM: Can't subset columns that don't exist
        test_id = test[id_cols_local].reset_index(drop=True)

        # Remove the outcome and ID columns from the training data
        train_x = train.drop([all_outcomes] + outcome_cols + id_cols_local, 
                             axis=1, errors="ignore")
        test_x = test.drop([all_outcomes] + outcome_cols + id_cols_local, 
                             axis=1, errors="ignore")

        # Set up and run the model
        
        if classifier == "lda":
            clf = sklearn.discriminant_analysis.LinearDiscriminantAnalysis()
        elif classifier == "nn":
            # Neural Net
            clf = sklearn.neural_network.MLPClassifier(random_state=1, 
                                                       max_iter=300)
        elif classifier == "svc":
            # Support Vector CLASSIFIER
            clf = sklearn.svm.SVC(kernel="linear", class_weight="balanced")
        else:
            raise ValueError(f"Unknown classifier requested: {classifier}")

        # Actually fit the model and do the prediction
        model = clf.fit(train_x, train_y)
        y_hat = model.predict(test_x)

        # Get the LD values and add it to the results
        if clf == "lda":
            y_transform = pd.DataFrame(model.transform(test_x))
            y_transform.columns = [f"LD{x}" for x 
                                    in range(1, y_transform.shape[1] + 1)]

            # print(f"Test data shape: {test_x.shape}")
            # print(f"Transform shape: {y_transform.shape}")

        # Put the predicted column back together with the ID cols
        results = pd.concat([test_id, test_y], axis=1)
        results["predicted"] = y_hat

        if clf == "lda":
            results = pd.concat([results, y_transform], axis=1)

        print(results)

        if not collapse_mcc:
            if not bootstrap_mcc:
                mcc1 = matthews_corrcoef(y_true=results[all_outcomes], 
                                        y_pred=results["predicted"])
                print(f"MCC:      {mcc1:.3f}")
                mcc_ci1 = [None, None]
            else:
                mcc1, mcc_ci1 = bootstrap_mcc_value(y_true=results[all_outcomes], 
                                                    y_pred=results["predicted"])

            # Create a vector of values
            group_row = pd.DataFrame({"group": hold_out, "mcc": float(mcc1),
                                      "mcc_ci_lo": float(mcc_ci1[0]),
                                      "mcc_ci_hi": float(mcc_ci1[1])},
                                     index=[0])

            group_mccs = pd.concat([group_mccs, group_row])

        # If requested, save the scalings to a CSV
        if dump_scalings and classifier == "lda":
            scalings = pd.DataFrame(model.scalings_)
            n_cols = scalings.shape[1]
            scalings.to_csv(scalings_name, index=False,
                            header=[f"LD{x}" for x in range(1, n_cols + 1)])

        results_df = pd.concat([results_df, results])

    # Show final results
    print(results_df)
    tab = pd.crosstab(results_df[all_outcomes], results_df["predicted"])

    # Display purposes only
    print("\nConfusion matrix:")
    print(tab)
    print()

    # Display basic accuracy
    acc = np.diag(tab).sum() / tab.to_numpy().sum()
    acc_cirange = (acc * (1 - acc) / tab.to_numpy().sum()) ** 0.5 * 1.96
    acc_ci = [acc - acc_cirange, acc + acc_cirange]
    print(f"Overall accuracy: {acc:.3f} [{acc_ci[0]:.3f}, {acc_ci[1]:.3f}]")

    if not bootstrap_mcc:
        # If not bootstrapping MCC, return missing values for the CI
        mcc = matthews_corrcoef(y_true=results_df[all_outcomes], 
                                y_pred=results_df["predicted"])
        print(f"MCC:      {mcc:.3f}")
        mcc_CI = [None, None]
    else:
        mcc, mcc_CI = bootstrap_mcc_value(y_true=results_df[all_outcomes], 
                                          y_pred=results_df["predicted"])

        print(f"Overall MCC: {mcc:.3f} [{mcc_CI[0]:.3f}, {mcc_CI[1]:.3f}]")

    # if len(outcome_cols) > 1:
    #     print(f"n\nStarting per-outcome calculations: {outcome_cols}")

    #     outcome = results_df["predicted"]
    #     outcome2 = outcome.str.split(".", expand=True)
    #     outcome2.columns = outcome_cols

    #     for o in outcome_cols:

    #         # print(f"\n{o}")
    #         print()
    #         tab = pd.crosstab(outcome2[o], data[o])
    #         print(tab)

    #         o_mcc, o_mcc_CI = bootstrap_mcc_value(y_true=data[o],
    #                                               y_pred=outcome2[o])

    #         print(f"MCC ({o}):   {o_mcc:.3f} [{o_mcc_CI[0]:.3f}, ",
    #               f"{o_mcc_CI[1]:.3f}]")

    results_df.to_csv(output_name, index=False)
    print(f"Saved to {output_name}!")

    # Put the overall MCC into an easily readable dataframe for later
    overall_mcc_df = pd.DataFrame({"mcc": float(mcc), 
                                   "mcc_ci_lo": float(mcc_CI[0]),
                                   "mcc_ci_hi": float(mcc_CI[1])},
                                  index=[0])

    return((results_df, group_mccs, overall_mcc_df))
