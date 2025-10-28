
import argparse as ap
import pickle
import pandas as pd
import os

parser = ap.ArgumentParser()

parser.add_argument("model")
parser.add_argument("data")
parser.add_argument("out")

args = parser.parse_args()

model_f = args.model
data_f = args.data
out = args.out

if not os.path.exists(out):

    with open(model_f, "rb") as f:
        print("Loading model")
        model = pickle.load(f)

    with open(data_f, "rb") as f:
        print("Loading data")
        data = pickle.load(f)

    print("Fitting")

    features = model.feature_names_in_
    data2 = data[features]
    id = data.drop(features, axis=1)

    xform = pd.DataFrame(model.transform(data2))
    predicted = pd.DataFrame(model.predict(data2))
    predicted.columns = ["predicted"]

    result = pd.concat([id, predicted, xform], axis=1)
    result.to_csv(out)

else:

    print(f"Output file {out} already exists, not running")