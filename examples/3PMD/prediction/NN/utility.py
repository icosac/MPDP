import numpy as np
import torch 
import pandas as pd
from sklearn.model_selection import train_test_split

def load_data(DATASET_PATH, using="cossin", val_prop=0.2, test_prop=0.2, scaler=None):
    data = pd.read_csv(DATASET_PATH, sep='\s+')

    features = ['kmax', 'theta_i', 'theta_f', 'alpha_m', 'alpha_f']
    target = 'th_m'

    if using == "cossin":
        # Transform angles into sine and cosine representations
        for col in ['theta_i', 'theta_f', 'alpha_m', 'alpha_f']:
            data[f'sin_{col}'] = np.sin(data[col])
            data[f'cos_{col}'] = np.cos(data[col])

        data['sin_th_m'] = np.sin(data['th_m'])
        data['cos_th_m'] = np.cos(data['th_m'])

        feature_columns = ['kmax'] + [f'sin_{col}' for col in ['theta_i', 'theta_f', 'alpha_m', 'alpha_f']] + [f'cos_{col}' for col in ['theta_i', 'theta_f', 'alpha_m', 'alpha_f']]
        X = data[feature_columns].values
        y = data[['sin_th_m', 'cos_th_m']].values
    elif using == "angle":
        X = data[features].values
        y = data[target].values
    elif using == "tan":
        for col in ['theta_i', 'theta_f', 'alpha_m', 'alpha_f']:
            data[f'tan_{col}'] = np.tan(data[col])

        data['tan_th_m'] = np.tan(data['th_m']/100.0)

        feature_columns = ['kmax'] + [f'tan_{col}' for col in ['theta_i', 'theta_f', 'alpha_m', 'alpha_f']]
        X = data[feature_columns].values
        y = data[['tan_th_m']].values


    X = scaler.fit_transform(X)

    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=val_prop, shuffle=True)
    if test_prop != 0:
        X_val, y_val, X_test, y_test, = train_test_split(X_val, y_val, test_size=test_prop, shuffle=True)

    X_train = torch.tensor(X_train, dtype=torch.float32)
    y_train = torch.tensor(y_train, dtype=torch.float32)
    X_val = torch.tensor(X_val, dtype=torch.float32)
    y_val = torch.tensor(y_val, dtype=torch.float32)
    if test_prop != 0:
        X_test = torch.tensor(X_test, dtype=torch.float32)
        y_test = torch.tensor(y_test, dtype=torch.float32)
    else:
        X_test = None
        y_test = None

    return X_train, y_train, X_val, y_val, X_test, y_test, scaler