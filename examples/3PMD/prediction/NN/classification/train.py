import torch
import copy
import numpy as np
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
import os
import matplotlib.pyplot as plt
from time import time as clock
from sklearn.metrics import precision_score, recall_score, accuracy_score, f1_score
from model import NeuralNet
from tqdm import tqdm  # Progress bar

from pathlib import Path

PROJECT_PATH = os.path.dirname(os.path.abspath(__file__))
DATASETS_PATH = os.path.join(Path(PROJECT_PATH).parent, "regression", "datasets")
DATASET_NAME = "small.csv"
DATASET_PATH = os.path.join(DATASETS_PATH, DATASET_NAME)
SAVE_PATH = os.path.join(PROJECT_PATH, "models")
SAVE_NAME = os.path.join(SAVE_PATH, "model.pt")
DO_PLOTS = True

# Set device to GPU if available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Using device:", device)

def model_train(model, X_train, y_train, X_val, y_val):
    model.to(device)

    # Loss function & optimizer
    loss_fn = nn.CrossEntropyLoss()  # For multi-class classification
    optimizer = optim.AdamW(model.parameters(), lr=0.001)

    n_epochs = 300
    batch_size = 8

    train_losses, val_losses = [], []
    train_accuracies, val_accuracies = [], []

    train_loader = DataLoader(TensorDataset(X_train, y_train), batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(TensorDataset(X_val, y_val), batch_size=batch_size, shuffle=False)

    best_acc = -np.inf
    best_weights = None

    for epoch in range(n_epochs):
        model.train()
        epoch_loss = 0.0
        correct = 0
        total = 0

        progress_bar = tqdm(train_loader, desc=f"Epoch {epoch+1}/{n_epochs}", dynamic_ncols=True)

        for X_batch, y_batch in progress_bar:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)

            # Convert one-hot to class indices
            y_batch = torch.argmax(y_batch, dim=1)

            optimizer.zero_grad()
            y_pred = model(X_batch)
            loss = loss_fn(y_pred, y_batch)
            loss.backward()
            optimizer.step()

            epoch_loss += loss.item() * X_batch.size(0)

            # Compute accuracy
            preds = torch.argmax(y_pred, dim=1)
            correct += (preds == y_batch).sum().item()
            total += y_batch.size(0)

            progress_bar.set_postfix(loss=loss.item())

        train_losses.append(epoch_loss / total)
        train_accuracies.append(correct / total)

        # Validation phase
        model.eval()
        val_loss = 0.0
        correct = 0
        total = 0

        with torch.no_grad():
            for X_batch, y_batch in val_loader:
                X_batch, y_batch = X_batch.to(device), y_batch.to(device)
                y_batch = torch.argmax(y_batch, dim=1)  # Convert one-hot to class indices

                y_pred = model(X_batch)
                loss = loss_fn(y_pred, y_batch)
                val_loss += loss.item() * X_batch.size(0)

                preds = torch.argmax(y_pred, dim=1)
                correct += (preds == y_batch).sum().item()
                total += y_batch.size(0)

        val_losses.append(val_loss / total)
        val_accuracies.append(correct / total)

        # Save best model
        if val_accuracies[-1] > best_acc:
            best_acc = val_accuracies[-1]
            best_weights = copy.deepcopy(model.state_dict())

    # Restore best model
    model.load_state_dict(best_weights)

    # Plotting
    if DO_PLOTS:
        plt.figure(figsize=(12, 6))
        
        plt.subplot(1, 2, 1)
        plt.plot(train_losses, label="Train Loss")
        plt.plot(val_losses, label="Validation Loss")
        plt.xlabel("Epoch")
        plt.ylabel("Loss")
        plt.legend()
        
        plt.subplot(1, 2, 2)
        plt.plot(train_accuracies, label="Train Accuracy")
        plt.plot(val_accuracies, label="Validation Accuracy")
        plt.xlabel("Epoch")
        plt.ylabel("Accuracy")
        plt.legend()
        
        plt.show()

    return best_acc


def evaluate_model(model, X_test, y_test):
    model.eval()
    with torch.no_grad():
        X_test = X_test.to(device)
        y_test = torch.argmax(y_test, dim=1).to(device)
        y_pred = torch.argmax(model(X_test), dim=1)
    return y_pred.cpu().numpy(), y_test.cpu().numpy()


def main():
    assert os.path.exists(DATASET_PATH), f"Dataset {DATASET_PATH} does not exist"
    if SAVE_NAME:
        os.makedirs(SAVE_PATH, exist_ok=True)

    time_start = clock()

    data = pd.read_csv(DATASET_PATH, sep='\s+')
    features = ['kmax', 'theta_i', 'theta_f', 'alpha_m', 'alpha_f']
    target = 'id_man_comb'

    X = data[features].values
    y = data[target].values

    # Encoding target labels
    encoder = OneHotEncoder(sparse_output=False)
    y = encoder.fit_transform(y.reshape(-1, 1))

    X = torch.tensor(X, dtype=torch.float32)
    y = torch.tensor(y, dtype=torch.float32)

    # Split data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y, random_state=42)

    # Train model
    model = NeuralNet()
    model_train(model, X_train, y_train, X_test, y_test)

    # Evaluate model
    y_pred_np, y_test_np = evaluate_model(model, X_test, y_test)

    # Print metrics
    print("Accuracy:", accuracy_score(y_test_np, y_pred_np))
    print("Precision:", precision_score(y_test_np, y_pred_np, average='weighted'))
    print("Recall:", recall_score(y_test_np, y_pred_np, average='weighted'))
    print("F1 Score:", f1_score(y_test_np, y_pred_np, average='weighted'))

    # Save model
    if SAVE_NAME:
        torch.save(model.state_dict(), SAVE_NAME)


if __name__ == "__main__":
    main()
