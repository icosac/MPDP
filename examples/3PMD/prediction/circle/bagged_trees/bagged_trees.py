from sklearn.ensemble import BaggingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score, f1_score
from sklearn.model_selection import train_test_split

import time, os
import pandas as pd
import numpy as np


print("Reading ds")
clock = time.time()

DATASETS_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "datasets")

data = pd.read_csv(os.path.join(DATASETS_PATH, 'small.csv'), sep='\s+')
# data = pd.read_csv(os.path.join(DATASETS_PATH, 'big_smaller.csv'), sep='\s+')
# data = pd.read_csv(os.path.join(DATASETS_PATH, 'big.csv'), sep='\s+')

features = ['kmax', 'theta_i', 'theta_f', 'alpha_m', 'alpha_f']
target = 'id_man_comb'

X = data[features].values
y = data[target].values

# Convert all X and y to float
X = X.astype(float)
y = y.astype(float)
print(f"Time to read ds: {time.time() - clock}")


print("Splitting ds")
clock = time.time()

reduce = False
if reduce:
    New_X, _, new_y, _ = train_test_split(X, y, test_size=1, random_state=42, shuffle=True)
else:
    New_X = X
    new_y = y

# Split the dataset into training (60%), validation (20%), and test (20%) sets
X_train, X_temp, y_train, y_temp = train_test_split(New_X, new_y, test_size=0.4, shuffle=True)
X_val, X_test, y_val, y_test = train_test_split(X_temp, y_temp, test_size=0.5, shuffle=True)

print(f"Time to split ds: {time.time() - clock}")

print("Training")
clock = time.time()

# Train a bagged tree ensemble
bagged_trees = BaggingClassifier(
    DecisionTreeClassifier(),
    n_estimators=51,
    n_jobs=8
)
bagged_trees.fit(X_train, y_train)
print(f"Time to train: {time.time() - clock}")


print("Predicting")
clock = time.time()
pred = bagged_trees.predict(X_test)
print(f"Time to predict: {time.time() - clock}")

print(f"Accuracy score: {accuracy_score(y_test, pred)}")
print(f"F1 score: {f1_score(y_test, pred, average='weighted')}")


# Exporting to C++
from skl2onnx import convert_sklearn
from skl2onnx.common.data_types import FloatTensorType

initial_type = [('float_input', FloatTensorType([None, X.shape[1]]))]
onnx_model = convert_sklearn(bagged_trees, initial_types=initial_type)

with open("bagged_trees.onnx", "wb") as f:
    f.write(onnx_model.SerializeToString())


