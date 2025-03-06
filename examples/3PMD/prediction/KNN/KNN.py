import pandas as pd
import numpy as np

import time

from annoy import AnnoyIndex

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score

print("Reading ds")
clock = time.time()
data = pd.read_csv('/home/enrico/Projects/mpdp/big.csv', sep='\s+')

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
    New_X, _, new_y, _ = train_test_split(X, y, test_size=0.999, shuffle=True)
else:
    New_X = X
    new_y = y

# Split the dataset into training (60%), validation (20%), and test (20%) sets
X_train, X_temp, y_train, y_temp = train_test_split(New_X, new_y, test_size=0.2, shuffle=True)
X_val, X_test, y_val, y_test = train_test_split(X_temp, y_temp, test_size=0.001, shuffle=True)

print(f"Time to split ds: {time.time() - clock}")

print(X_train.shape, X_val.shape, X_test.shape)


training = True
annoy_index = AnnoyIndex(5, 'euclidean')

if training:
    print("Training")
    clock = time.time()
    for i in range(len(X_train)):
        annoy_index.add_item(i, X_train[i])
    annoy_index.build(11, n_jobs=8)
    print(f"Time to train: {time.time() - clock}")

    print("Saving index")
    clock = time.time()
    annoy_index.save('knn.ann')
    print(f"Time to save index: {time.time() - clock}")

    # Save X_val and y_val
    np.save('X_val.npy', X_val)
    np.save('y_val.npy', y_val)
    np.save('y_train.npy', y_train)

    
else:
    print("Loading index")
    clock = time.time()
    annoy_index.load('knn.ann')
    print(f"Time to load index: {time.time() - clock}")

    # Load X_val and y_val
    X_val = np.load('X_val.npy')
    y_val = np.load('y_val.npy')
    y_train = np.load('y_train.npy')


def predict_knn(X_val, k=3):
    """Predict labels for X_val using k-nearest neighbors with majority voting."""
    y_pred = []
    
    for x in X_val:
        indices = annoy_index.get_nns_by_vector(x, k)  # Get k nearest neighbors
        nearest_labels = y_train[indices].astype(int)  # Retrieve corresponding labels
        predicted_label = np.bincount(nearest_labels).argmax()  # Majority vote
        y_pred.append(predicted_label)
    
    return np.array(y_pred)


print("Predicting")
clock = time.time()

# Make predictions
y_pred = predict_knn(X_val, k=5)

# Evaluate performance
accuracy = accuracy_score(y_val, y_pred)
f1 = f1_score(y_val, y_pred, average='weighted')

print(f"Accuracy: {accuracy*100:.1f}%")
print(f"F1 Score: {f1:.4f}")

print(f"Time to predict: {time.time() - clock}")








