import torch
import copy
import numpy as np
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import pandas as pd
from sklearn.model_selection import train_test_split
import os
import matplotlib.pyplot as plt
from time import time as clock
from model import NeuralNet
from tqdm import tqdm
from pathlib import Path
from sklearn.preprocessing import StandardScaler

PROJECT_PATH = os.path.dirname(os.path.abspath(__file__))
DATASETS_PATH = os.path.join(Path(PROJECT_PATH).parent.parent, "datasets")
DATASET_NAME = "small.csv"
DATASET_PATH = os.path.join(DATASETS_PATH, DATASET_NAME)
SAVE_PATH = os.path.join(PROJECT_PATH, "models")
SAVE_NAME = os.path.join(SAVE_PATH, "model1.pt")

DO_PLOTS = True

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# device = "cpu"
print("Using device:", device)

def angular_loss(y_pred, y_true):
    """Custom loss function for periodic angles."""
    sin_loss = nn.MSELoss()(y_pred[:, 0], y_true[:, 0])
    cos_loss = nn.MSELoss()(y_pred[:, 1], y_true[:, 1])
    return sin_loss + cos_loss

def model_train(model, X_train, y_train, X_val, y_val):
    model.to(device)
    
    optimizer = optim.Adam(model.parameters(), lr=0.001)  # Increased learning rate
    n_epochs = 10  
    batch_size = 32

    train_losses, val_losses = [], []

    train_loader = DataLoader(TensorDataset(X_train, y_train), batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(TensorDataset(X_val, y_val), batch_size=batch_size, shuffle=False)

    best_loss = float("inf")
    best_weights = None

    for epoch in range(n_epochs):
        model.train()
        epoch_loss = 0.0
        progress_bar = tqdm(train_loader, desc=f"Epoch {epoch+1}/{n_epochs}", dynamic_ncols=True)

        for X_batch, y_batch in progress_bar:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)

            optimizer.zero_grad()
            y_pred = model(X_batch)
            loss = angular_loss(y_pred, y_batch)

            loss.backward()
            optimizer.step()

            epoch_loss += loss.item() * len(X_batch)
            progress_bar.set_postfix(loss=loss.item())

        train_loss = epoch_loss / len(X_train)
        train_losses.append(train_loss)

        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for X_batch, y_batch in val_loader:
                X_batch, y_batch = X_batch.to(device), y_batch.to(device)
                y_pred = model(X_batch)
                loss = angular_loss(y_pred, y_batch)
                val_loss += loss.item() * len(X_batch)

        val_loss /= len(X_val)
        val_losses.append(val_loss)

        if val_loss < best_loss:
            best_loss = val_loss
            best_weights = copy.deepcopy(model.state_dict())

        # Print accuracy and MAE every 2 epochs
        if (epoch + 1) % 2 == 0:
            y_val_pred = evaluate_model(model, X_val)
            actual_angles = np.arctan2(y_val[:, 0].cpu().numpy(), y_val[:, 1].cpu().numpy())
            error = np.abs(actual_angles - y_val_pred)
            mae = np.mean(error)  # Mean Absolute Error
            print(f"Epoch {epoch+1}: Mean Absolute Error: {mae:.4f}")

            save_path_epoch = os.path.join(SAVE_PATH, f"model_epoch_{epoch+1}.pt")
            torch.save(model.state_dict(), save_path_epoch)
            print(f"Model saved at {save_path_epoch}")

    model.load_state_dict(best_weights)

    if DO_PLOTS:
        plt.figure(figsize=(10, 5))
        plt.plot(train_losses, label="Train Loss")
        plt.plot(val_losses, label="Validation Loss")
        plt.xlabel("Epoch")
        plt.ylabel("Loss")
        plt.title("Training vs Validation Loss")
        plt.legend()
        # plt.show()
        plt.savefig(f"regression{n_epochs}-{batch_size}.png")

    return best_loss

def evaluate_model(model, X_test):
    model.eval()
    with torch.no_grad():
        X_test = X_test.to(device)
        y_pred = model(X_test).cpu().numpy()
    return np.arctan2(y_pred[:, 0], y_pred[:, 1])  # Convert (sin, cos) back to angle

def main():
    assert os.path.exists(DATASET_PATH), f"Dataset {DATASET_PATH} does not exist"
    os.makedirs(SAVE_PATH, exist_ok=True)

    time_start = clock()
    
    data = pd.read_csv(DATASET_PATH, sep='\s+')

    features = ['kmax', 'theta_i', 'theta_f', 'alpha_m', 'alpha_f']
    target = 'th_m'

    # Transform angles into sine and cosine representations
    for col in ['theta_i', 'theta_f', 'alpha_m', 'alpha_f']:
        data[f'sin_{col}'] = np.sin(data[col])
        data[f'cos_{col}'] = np.cos(data[col])

    data['sin_th_m'] = np.sin(data['th_m'])
    data['cos_th_m'] = np.cos(data['th_m'])

    feature_columns = ['kmax'] + [f'sin_{col}' for col in ['theta_i', 'theta_f', 'alpha_m', 'alpha_f']] + [f'cos_{col}' for col in ['theta_i', 'theta_f', 'alpha_m', 'alpha_f']]
    X = data[feature_columns].values
    y = data[['sin_th_m', 'cos_th_m']].values

    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, shuffle=True)

    X_train = torch.tensor(X_train, dtype=torch.float32)
    y_train = torch.tensor(y_train, dtype=torch.float32)
    X_test = torch.tensor(X_test, dtype=torch.float32)
    y_test = torch.tensor(y_test, dtype=torch.float32)

    model = NeuralNet(in_size=len(feature_columns))
    model_train(model, X_train, y_train, X_test, y_test)

    y_pred = evaluate_model(model, X_test)

    actual_angles = np.arctan2(y_test[:, 0], y_test[:, 1])  # Convert ground truth to angle
    error = np.abs(actual_angles - y_pred)
    # error = np.minimum(error, 2 * np.pi - error)  # Handle circular distance

    print("Mean Absolute Error (Angular):", np.average(error))

    if SAVE_NAME:
        torch.save(model.state_dict(), SAVE_NAME)

if __name__ == "__main__":
    main()
