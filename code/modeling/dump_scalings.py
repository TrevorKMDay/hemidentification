
import argparse as ap
import pickle
import pandas as pd
import os

parser = ap.ArgumentParser()

parser.add_argument("pickle",
                    help="Column to predict")

parser.add_argument("out")

args = parser.parse_args()

file = args.pickle
out = args.out

if not os.path.exists(out):

    with open(file, "rb") as f:
        data = pickle.load(f)

    features = pd.DataFrame(data.feature_names_in_)
    features.columns = ["feature"]

    scalings = pd.DataFrame(data.scalings_)
    scalings.columns = ["LD" + str(x + 1) for x in range(len(scalings.columns))]

    x = pd.concat([features.reset_index(drop=True),
                   scalings.reset_index(drop=True)],
                  axis=1)

    x.to_csv(out, index=False)