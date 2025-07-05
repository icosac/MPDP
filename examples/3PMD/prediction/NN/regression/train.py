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
from tqdm import tqdm  # Import tqdm for progress bar
from pathlib import Path
from sklearn.metrics import mean_squared_error, mean_absolute_error

PROJECT_PATH = os.path.dirname(os.path.abspath(__file__))

DATASETS_PATH = os.path.join(Path(PROJECT_PATH).parent, "..", "datasets")
DATASET_NAME = "small.csv"
DATASET_PATH = os.path.join(DATASETS_PATH, DATASET_NAME)

SAVE_PATH = os.path.join(PROJECT_PATH, "models")
SAVE_NAME = os.path.join(SAVE_PATH, "model.pt")

DO_PLOTS = True

# Set device to GPU if available, otherwise use CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# device = "cpu"
print("Using device:", device)

def model_train(model, X_train, y_train, X_val, y_val):
    model.to(device)

    # Use MSE loss for regression
    loss_fn = nn.MSELoss()  
    optimizer = optim.Adam(model.parameters(), lr=0.0001)

    n_epochs = 5
    batch_size = 16  

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
            loss = loss_fn(y_pred, y_batch)

            loss.backward()
            optimizer.step()

            epoch_loss += loss.item() * len(X_batch)
            progress_bar.set_postfix(mse=loss.item())

        train_loss = epoch_loss / len(X_train)
        train_losses.append(train_loss)

        # Validation
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for X_batch, y_batch in val_loader:
                X_batch, y_batch = X_batch.to(device), y_batch.to(device)
                y_pred = model(X_batch)
                loss = loss_fn(y_pred, y_batch)
                val_loss += loss.item() * len(X_batch)

        val_loss /= len(X_val)
        val_losses.append(val_loss)

        if val_loss < best_loss:
            best_loss = val_loss
            best_weights = copy.deepcopy(model.state_dict())

    # Restore the best model
    model.load_state_dict(best_weights)

    if DO_PLOTS:
        plt.figure(figsize=(10, 5))
        plt.plot(train_losses, label="Train MSE")
        plt.plot(val_losses, label="Validation MSE")
        plt.xlabel("Epoch")
        plt.ylabel("Loss (MSE)")
        plt.title("Training vs Validation Loss")
        plt.legend()
        # plt.show()
        plt.savefig("loss_regression.png")

    return best_loss

# Function to evaluate the model
def evaluate_model(model, X_test, y_test):
    model.eval()
    with torch.no_grad():
        X_test = X_test.to(device)
        y_pred = model(X_test)
    return y_pred.cpu().numpy().flatten()

def main():
    assert os.path.exists(DATASET_PATH), f"Dataset {DATASET_PATH} does not exist"
    os.makedirs(SAVE_PATH, exist_ok=True)

    time_start = clock()
    
    data = pd.read_csv(DATASET_PATH, sep='\s+')

    features = ['kmax', 'theta_i', 'theta_f', 'alpha_m', 'alpha_f']
    target = 'th_m'

    X = data[features].values
    y = data[target].values

    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, shuffle=True)

    X_train = torch.tensor(X_train, dtype=torch.float32)
    y_train = torch.tensor(y_train, dtype=torch.float32).reshape(-1, 1)
    X_test = torch.tensor(X_test, dtype=torch.float32)
    y_test = torch.tensor(y_test, dtype=torch.float32).reshape(-1, 1)

    # Train the model
    model = NeuralNet()
    model_train(model, X_train, y_train, X_test, y_test)

    # Evaluate the model
    y_pred = evaluate_model(model, X_test, y_test.cpu().numpy())

    print("Mean Squared Error:", mean_squared_error(y_test.cpu().numpy(), y_pred))
    print("Mean Absolute Error:", mean_absolute_error(y_test.cpu().numpy(), y_pred))

    if SAVE_NAME:
        torch.save(model.state_dict(), SAVE_NAME)

if __name__ == "__main__":
    main()
