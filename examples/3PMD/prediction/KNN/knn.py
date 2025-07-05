import pandas as pd
import numpy as np

import time, os

from annoy import AnnoyIndex

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score


RUN_TESTS       = 1 # Set 1 or 2 or 3 to run test set 1, 2, or both
DATASETS_PATH   = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "datasets")
FEATURES        = ['kmax', 'theta_i', 'theta_f', 'alpha_m', 'alpha_f']
TARGET          = 'id_man_comb'
TRAINING        = True

annoy_index = AnnoyIndex(5, 'euclidean')

if TRAINING:
    print("###############\nTraining")
    print("Reading ds")
    clock = time.time()

    # data = pd.read_csv(os.path.join(DATASETS_PATH, 'small.csv'), sep='\s+')
    data = pd.read_csv(os.path.join(DATASETS_PATH, 'big_smaller_new.csv'), sep='\s+')
    # data = pd.read_csv(os.path.join(DATASETS_PATH, 'big.csv'), sep='\s+')
    
    # Order data for the id_man_comb column so that the labels are in increasing order
    data = data.sort_values("id_man_comb", axis=0, ascending=True)

    X = data[FEATURES].values
    y = data[TARGET].values

    # Convert all X and y to float
    # X = X.astype(float)
    # y = y.astype(float)
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
    # X_train, X_temp, y_train, y_temp = train_test_split(New_X, new_y, test_size=0.1, shuffle=True)
    # if X_temp is not None:
    #     X_val, X_test, y_val, y_test = train_test_split(X_temp, y_temp, test_size=0.2, shuffle=True)
    # else:
    X_val = None
    X_test = None
    y_val = None
    y_test = None

    X_train = X
    y_train = y

    print(f"Time to split ds: {time.time() - clock}")    
    
    clock = time.time()
    for i in range(len(X_train)):
        annoy_index.add_item(i, X_train[i])
    annoy_index.build(11, n_jobs=8)
    print(f"Time to train: {time.time() - clock}")

    print("Saving index")
    clock = time.time()
    annoy_index.save('knn_new.ann')
    print(f"Time to save index: {time.time() - clock}")

    # Save X_val and y_val
    # np.savetxt('X_val_new.csv', np.asarray(X_val), delimiter=' ')
    # np.savetxt('y_val_new.csv', np.asarray(y_val), delimiter=' ')
    np.savetxt('y_train_new.csv', np.asarray(y_train), delimiter=' ')
    
    map_labels = {y_train[0]: [0]}
    label = y_train[0]
    for i in range(len(y_train)-1):
        if label != y_train[i+1]:
            map_labels[label].append(i)
            map_labels[y_train[i+1]] = [i+1]
            label = y_train[i+1]
    map_labels[label].append(len(y_train)-1)
    with open("y_labels_intervals_new.csv", "w") as f:
        for k, v in map_labels.items():
            f.write(f"{k} {v[0]} {v[1]}\n")
    

    
else:
    print("###############\nLoading index")
    clock = time.time()
    annoy_index.load('knn.ann')
    print(f"Time to load index: {time.time() - clock}")

    # Load X_val and y_val
    # X_val = np.load('X_val.npy')
    # y_val = np.load('y_val.npy')
    # y_train = np.load('y_train.npy')

    # X_val = np.loadtxt('X_val.csv', delimiter=' ')
    # y_val = np.loadtxt('y_val.csv', delimiter=' ')
    y_train = np.loadtxt('y_train.csv', delimiter=' ')


def predict_knn(X_val, k=3):
    """Predict labels for X_val using k-nearest neighbors with majority voting."""
    y_pred = []
    
    for x in X_val:
        indices = annoy_index.get_nns_by_vector(x, k)  # Get k nearest neighbors
        nearest_labels = y_train[indices].astype(int)  # Retrieve corresponding labels
        # for i, l in zip (indices, nearest_labels):
        #     print(f"Index: {i}, Label: {l}")
        predicted_label = np.bincount(nearest_labels).argmax()  # Majority vote
        y_pred.append(predicted_label)
    
    return np.array(y_pred)


print("###############\nPredicting")
# y_pred = predict_knn([[0.5, 0.2, 0.8, 0.3, 0.9]], k=4)
# print(y_pred)
# import sys
# sys.exit(0)

if X_val is not None:
    print("###############\nPredicting validation set")
    clock = time.time()
    # Make predictions
    y_pred = predict_knn(X_val, k=4)

    # Evaluate performance
    accuracy = accuracy_score(y_val, y_pred)
    f1 = f1_score(y_val, y_pred, average='weighted')

    print(f"Accuracy: {accuracy*100:.1f}%")
    print(f"F1 Score: {f1:.4f}")

    print(f"Time to predict: {time.time() - clock}")

if X_test is not None:
    print("###############\nPredicting test set")
    clock = time.time()
    # Make predictions
    y_pred = predict_knn(X_test, k=4)

    # Evaluate performance
    accuracy = accuracy_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred, average='weighted')

    print(f"Accuracy: {accuracy*100:.1f}%")
    print(f"F1 Score: {f1:.4f}")

    print(f"Time to predict: {time.time() - clock}")

if RUN_TESTS in [1, 3]:
    print("###############\nRunning test1")

    # data = pd.read_csv(os.path.join(DATASETS_PATH, 'testset1.csv'), sep='\s+')
    data = pd.read_csv("/Users/enrico/Projects/mpdp/ds_out_4.csv", sep='\s+')
    X_test = data[FEATURES].values
    y_test = data[TARGET].values

    clock = time.time()
    # Make predictions
    y_pred = predict_knn(X_test, k=4)

    # Evaluate performance
    accuracy = accuracy_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred, average='weighted')

    print(f"Accuracy: {accuracy*100:.1f}%")
    print(f"F1 Score: {f1:.4f}")
    print(f"Time to predict: {time.time() - clock}")

if RUN_TESTS in [2, 3]:
    print("###############\nRunning test2")

    data = pd.read_csv(os.path.join(DATASETS_PATH, 'testset2.csv'), sep='\s+')
    X_test = data[FEATURES].values
    y_test = data[TARGET].values

    clock = time.time()
    # Make predictions
    y_pred = predict_knn(X_test, k=4)

    # Evaluate performance
    accuracy = accuracy_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred, average='weighted')

    print(f"Accuracy: {accuracy*100:.1f}%")
    print(f"F1 Score: {f1:.4f}")
    print(f"Time to predict: {time.time() - clock}")