import pickle
import pandas as pd

from pathlib import Path
from imblearn.over_sampling import RandomOverSampler

# Load data

home = Path("/Users/tkmd/Google Drive/My Drive/Projects/"
            "hemisphere_fingerprinting/code/modeling/inputs/")

with open(home / "hemiconnectome.pickle", "rb") as f:
    data = pickle.load(f)

print(f"Original shape: {data.shape}")

# Create oversampling

ros = RandomOverSampler(random_state=0,
                        sampling_strategy="minority")

data_y = data["handedness"]
data_x = data.drop(["handedness"], axis=1)

new_x, new_y = ros.fit_resample(data_x, data_y)

print(f"New shape: X: {new_x.shape}; Y: {new_y.shape}")
# print(new_y.value_counts())

new_data = pd.concat([new_y, new_x], axis=1)

# print(new_data.handedness.value_counts())
# print(new_data)
# print(new_data["sub"].value_counts())

new_data.to_pickle("hemiconnectome_oversampled.pickle")
