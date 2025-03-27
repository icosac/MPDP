import torch
import copy
import numpy as np
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import os, sys
import matplotlib.pyplot as plt
from model import NeuralNet
from tqdm import tqdm
from pathlib import Path
from sklearn.preprocessing import StandardScaler
import joblib  # Add this import for saving the scaler
import json  # Add this import for saving scaler parameters as JSON

# import utility which is one folder back
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from utility import load_data

PROJECT_PATH = os.path.dirname(os.path.abspath(__file__))
DATASETS_PATH = os.path.join(Path(PROJECT_PATH).parent.parent, "datasets")
DATASET_NAME = "small.csv"
DATASET_PATH = os.path.join(DATASETS_PATH, DATASET_NAME)
SAVE_PATH = os.path.join(PROJECT_PATH, "models")
SAVE_NAME = os.path.join(SAVE_PATH, "model1.pt")

DO_PLOTS = True

USING = "tan"  # "cossin", "angle", "tan"

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# device = "cpu"
print("Using device:", device)

def loss_fn(y_pred, y_true):
    """Custom loss function for angles."""
    if USING == "cossin":
        return cossin_loss(y_pred, y_true)
    elif USING == "angle":
        return angle_loss(y_pred, y_true)
    elif USING == "tan":
        return tan_loss(y_pred, y_true)

def cossin_loss(y_pred, y_true):
    """Custom loss function for periodic angles."""
    sin_loss = nn.MSELoss()(y_pred[:, 0], y_true[:, 0])
    cos_loss = nn.MSELoss()(y_pred[:, 1], y_true[:, 1])
    return sin_loss + cos_loss

def angle_loss(y_pred, y_true):
    """Custom loss function for angles."""
    return nn.MSELoss()(y_pred, y_true)

def tan_loss(y_pred, y_true):
    """Custom loss function for angles."""
    return nn.MSELoss()(y_pred, y_true)

def convert_labels(y):
    if USING == "cossin":
        actual_labels = np.arctan2(y[:, 0], y[:, 1]) # Convert ground truth to angle
    elif USING == "angle":
        actual_labels = y
    elif USING == "tan":
        actual_labels = np.arctan(y) # Convert ground truth to angle
    return actual_labels.reshape(-1)

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
            loss = loss_fn(y_pred, y_batch)

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
                loss = loss_fn(y_pred, y_batch)
                val_loss += loss.item() * len(X_batch)

        val_loss /= len(X_val)
        val_losses.append(val_loss)

        if val_loss < best_loss:
            best_loss = val_loss
            best_weights = copy.deepcopy(model.state_dict())

        # Print accuracy and MAE every 2 epochs
        if (epoch + 1) % 2 == 0:
            y_val_pred = evaluate_model(model, X_val)
            actual_labels = convert_labels(y_val)
            error = np.abs(actual_labels - y_val_pred).cpu().numpy()
            error = np.minimum(error, 2 * np.pi - error)  # Handle circular distance
            mae = np.mean(error)  # Mean Absolute Error
            print(f"Epoch {epoch+1}: Mean Absolute Error: {mae}")

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
    return convert_labels(y_pred) 


def main():
    assert os.path.exists(DATASET_PATH), f"Dataset {DATASET_PATH} does not exist"
    os.makedirs(SAVE_PATH, exist_ok=True)
    
    scaler = StandardScaler()
    X_train, y_train, X_val, y_val, _, _, scaler = load_data(DATASET_PATH, using=USING, scaler=scaler, val_prop=0.2, test_prop=0.0)

    n_features = X_train.shape[1]
    n_labels = 1 if len(y_val.shape) == 1 else y_val.shape[1]

    model = NeuralNet(in_size=n_features, out_size=n_labels)
    model_train(model, X_train, y_train, X_val, y_val)

    y_pred = evaluate_model(model, X_val)

    actual_labels = convert_labels(y_val)

    error = np.abs(actual_labels - y_pred)
    error = np.minimum(error, 2 * np.pi - error)  # Handle circular distance  

    print("Mean Absolute Error (Angular):", np.average(error))

    if SAVE_NAME:
        # Saving model for Python3
        torch.save(model.state_dict(), SAVE_NAME)
        # Saving model to torchscript for C++
        model.to(device)
        model.eval()
        example = torch.rand(1, n_features).to(device)
        traced_script_module = torch.jit.trace(model, example)
        traced_script_module.save(SAVE_NAME.replace(".pt", "_cpp.pt"))

        # Save the scaler
        scaler_save_path = SAVE_NAME.replace(".pt", "_scaler.pkl")
        joblib.dump(scaler, scaler_save_path)
        print(f"Scaler saved at {scaler_save_path}")

        # Save scaler parameters as JSON
        scaler_params = {"mean": scaler.mean_.tolist(), "scale": scaler.scale_.tolist()}
        scaler_json_path = SAVE_NAME.replace(".pt", "_scaler.json")
        with open(scaler_json_path, "w") as f:
            json.dump(scaler_params, f)
        print(f"Scaler parameters saved at {scaler_json_path}")

if __name__ == "__main__":
    main()
