#/usr/bin/env python3

from math import pi, cos, sin, radians, degrees
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 

kmax = 1

theta_i = radians(0)
theta_f = radians(180)
alpha_m = radians(45)
alpha_f = radians(90)

pi_c = (1, 0, theta_i)
pm_c = (cos(alpha_m), sin(alpha_m), 0)
pf_c = (cos(alpha_f), sin(alpha_f), theta_f)

SEP = '|'

FULL_PATH_DIR = os.path.dirname(os.path.realpath(__file__))
COMBINATIONS_FILE = os.path.join(FULL_PATH_DIR, 'combinations.csv')
RESULTS_FILE = os.path.join(FULL_PATH_DIR, 'results.txt')
LOG_FILE = os.path.join(FULL_PATH_DIR, 'log.txt')

TENSORFLOW = True

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

data = pd.read_csv('/home/enrico/Projects/mpdp/3PDS_Circle32_1_1_1.csv', sep='\s+')

features = ['kmax', 'theta_i', 'theta_f', 'alpha_m', 'alpha_f']
target = ['th_m']

X = data[features].values
y = data[target].values

# Split the dataset into training (60%), validation (20%), and test (20%) sets
X_train, X_temp, y_train, y_temp = train_test_split(X, y, test_size=0.4, random_state=42, shuffle=True)
X_val, X_test, y_val, y_test = train_test_split(X_temp, y_temp, test_size=0.5, random_state=42, shuffle=True)

# Standardize the features
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_val = scaler.transform(X_val)
X_test = scaler.transform(X_test)

import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.optimizers import Adam

print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))

from datetime import datetime
import itertools

def predict_id_man_comb(model, kmax, theta_i, theta_f, alpha_m, alpha_f, scaler = None):
    input_data = np.array([[kmax, theta_i, theta_f, alpha_m, alpha_f]])
    if scaler:
        input_data = scaler.transform(input_data)

    [th] = model.predict(input_data)
    
    print(th)
    return th

def read_remaining_combinations():
    try:
        data = pd.read_csv(COMBINATIONS_FILE, sep=SEP, header=None, engine='python')
    except pd.errors.EmptyDataError:
        print("No more combinations")
        return None
    
    failed = False

    try:
        tmp_values = data.values[0][0].split('(')[1]
        tmp_values = tmp_values.split(')')[0]
        tmp_values = tmp_values.split(', ')
        next_combination = [int(x) for x in tmp_values]
    except:
        failed = True
        print(f"Skipping {data.values[0][0]}")
        
    if failed:
        return None

    remaining_combinations = list(data.values[0][1:])
    print(f"Remaining combinations: {len(remaining_combinations)}")
    write_remaining_combinations(remaining_combinations)

    return next_combination 


def write_remaining_combinations(combinations : list):
    with open(COMBINATIONS_FILE, 'w') as f:
        for comb_id in range(len(combinations)):
            if comb_id > 0:
                f.write(f'{SEP}')
            if combinations[comb_id] is not None:
                f.write(f'{combinations[comb_id]}')
                

def first_run():
    if not os.path.exists(RESULTS_FILE):
        with open(RESULTS_FILE, 'w') as f:
            f.write('n_layers,test_loss,pred_th,neurons\n')

    if os.path.exists(COMBINATIONS_FILE):
        res = input('File already exists, do you want to overwrite the file? (y/N) ')
        if res.lower() != 'y':
            print("Quitting")
            return

    num_layers = range(4, 6)
    num_neurons = [32, 64, 128, 1024, 4096]

    combinations = []
    for num_layer in num_layers:
        combinations += list(itertools.product(num_neurons, repeat=num_layer))

    print(f"Testing {len(combinations)} combinations")

    write_remaining_combinations(combinations)
    

def train_and_predict(combination):
    n_layers = len(combination)
    neurons = combination

    print(f'layers: {n_layers}, neurons: {neurons}')
    model = Sequential()
    model.add(Dense(neurons[0], activation='relu', input_shape=(X_train.shape[1],)))
    for n_neurons in neurons[1:-1]:
        model.add(Dense(n_neurons, activation='relu'))
    model.add(Dense(1))

    model.compile(optimizer=Adam(learning_rate=0.01), loss='mse')

    t = datetime.now()
    model.fit(X_train, y_train, validation_data=(X_val, y_val), epochs=50, batch_size=128, verbose=1)
    train_time = (datetime.now() - t)/1000.0

    test_loss = model.evaluate(X_test, y_test)
    # print(f'Test Mean Squared Error: {test_loss}') 

    t = datetime.now()
    pred = predict_id_man_comb(model, kmax, theta_i, theta_f, alpha_m, alpha_f, scaler)
    pred_time = (datetime.now() - t)/1000.0
    # print(f"elapsed {pred_time}")
    with open(RESULTS_FILE, 'a') as f:
        if pred is not None:
        # print(f'The angle of the middle point is: {degrees(pred)}, {pred}')
        # print(f'The angle of the middle point is: {THREE_P_MAN_ID[int(man)-1]}')
            f.write(f'{n_layers},{test_loss},{pred[0]},{pred_time},{train_time},{neurons}\n')
        else:
            f.write(f'{n_layers},{test_loss},None,None,{pred_time},{train_time},{neurons}\n')




def main():
    # train_and_predict()
    if False:
        first_run()
    else:
        combination = read_remaining_combinations()
        print(f'Testing {combination}')
        with open(LOG_FILE, "a") as f:
            f.write(f'Testing {combination}\n')

        train_and_predict(combination)
        os.system("sudo reboot")

if __name__ == '__main__':
    main()
