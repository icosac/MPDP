from sklearn.metrics import precision_score, recall_score, accuracy_score, f1_score, confusion_matrix, ConfusionMatrixDisplay
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
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from model import NeuralNet
from tqdm import tqdm  # Import tqdm for progress bar

from pathlib import Path

PROJECT_PATH = os.path.dirname(os.path.abspath(__file__))

DATASETS_PATH = os.path.join(Path(PROJECT_PATH).parent, "datasets")
DATASET_NAME = "small.csv"
DATASET_PATH = os.path.join(DATASETS_PATH, DATASET_NAME)

SAVE_PATH = os.path.join(PROJECT_PATH, "models")
SAVE_NAME = os.path.join(SAVE_PATH, "model.pt")

DO_PLOTS = True

# Set device to GPU if available, otherwise use CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Using device:", device)

def model_train(model, X_train, y_train, X_val, y_val):
    # Move model to the device (GPU or CPU)
    model.to(device)

    # Loss function and optimizer
    loss_fn = nn.CrossEntropyLoss()  # Binary Cross-Entropy Loss
    optimizer = optim.AdamW(model.parameters(), lr=0.0001)

    n_epochs = 40  # Number of epochs to run
    batch_size = 16  # Size of each batch

    # Metrics tracking
    train_losses, val_losses = [], []
    train_accuracies, val_accuracies = [], []

    train_loader = DataLoader(TensorDataset(X_train, y_train), batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(TensorDataset(X_val, y_val), batch_size=batch_size, shuffle=False)

    # Hold the best model
    best_acc = -np.inf
    best_weights = None

    for epoch in range(n_epochs):
        # Training loop with tqdm progress bar
        model.train()
        epoch_loss = 0.0
        correct = 0
        total = 0

        # Create the progress bar
        progress_bar = tqdm(train_loader, desc=f"Epoch {epoch+1}/{n_epochs}", dynamic_ncols=True)

        for X_batch, y_batch in progress_bar:
            # Move data to the device
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)

            # Forward pass
            y_pred = model(X_batch)
            loss = loss_fn(y_pred, y_batch)

            # Backward pass
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            # Accumulate metrics
            epoch_loss += loss.item() * len(X_batch)
            correct += (y_pred.round() == y_batch).float().sum().item()
            total += len(y_batch)

            # Update the progress bar with stats
            progress_bar.set_postfix(loss=epoch_loss / total, accuracy=correct / total)

        # Record training metrics
        train_loss = epoch_loss / total
        train_acc = correct / total
        train_losses.append(train_loss)
        train_accuracies.append(train_acc)

        # Validation loop
        model.eval()
        epoch_loss = 0.0
        correct = 0
        total = 0
        with torch.no_grad():
            for X_batch, y_batch in val_loader:
                # Move data to the device
                X_batch, y_batch = X_batch.to(device), y_batch.to(device)
                y_pred = model(X_batch)
                loss = loss_fn(y_pred, y_batch)
                epoch_loss += loss.item() * len(X_batch)
                correct += (y_pred.round() == y_batch).float().sum().item()
                total += len(y_batch)

        # Record validation metrics
        val_loss = epoch_loss / total
        val_acc = correct / total
        val_losses.append(val_loss)
        val_accuracies.append(val_acc)

        # Save the best model
        if val_acc > best_acc:
            best_acc = val_acc
            best_weights = copy.deepcopy(model.state_dict())

    # Restore the best model
    model.load_state_dict(best_weights)

    if DO_PLOTS:
        # Plot metrics
        plt.figure(figsize=(12, 6))

        # Plot loss
        plt.subplot(1, 2, 1)
        plt.plot(train_losses, label="Train Loss")
        plt.plot(val_losses, label="Validation Loss")
        plt.xlabel("Epoch")
        plt.ylabel("Loss")
        plt.title("Loss per Epoch")
        plt.legend()

        # Plot accuracy
        plt.subplot(1, 2, 2)
        plt.plot(train_accuracies, label="Train Accuracy")
        plt.plot(val_accuracies, label="Validation Accuracy")
        plt.xlabel("Epoch")
        plt.ylabel("Accuracy")
        plt.title("Accuracy per Epoch")
        plt.legend()

        plt.tight_layout()
        plt.show()

    return best_acc


# Function to evaluate the model
def evaluate_model(model, X_test, y_test):
    model.eval()  
    with torch.no_grad():
        X_test = X_test.to(device)  # Move test data to GPU
        y_pred = model(X_test).round()
    return y_pred

def main():
    assert(os.path.exists(DATASET_PATH)), f"Dataset {DATASET_PATH} doest not exist"
    if SAVE_NAME != "":
        os.makedirs(SAVE_PATH, exist_ok=True)

    time_start = clock()
    
    data = pd.read_csv(DATASET_PATH, sep='\s+')
    X = data.iloc[:, :6]    # features   
    y = data.iloc[:, 6]      # target

    # Encoding the target
    encoder = OneHotEncoder(sparse_output=False)  # Use dense array output
    y = encoder.fit_transform(y.to_numpy().reshape(-1, 1))  # Convert y to numpy and reshape

    X = torch.tensor(X.values, dtype=torch.float32)
    y = torch.tensor(y, dtype=torch.float32)

    # Train-test split: Hold out the test set for final model evaluation
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y, random_state=42)

    # Train the model
    model = NeuralNet()
    model_train(model, X_train, y_train, X_test, y_test)

    # Evaluate the model on the test set
    y_pred = evaluate_model(model, X_test, y_test)

    # Convert predictions and true labels to numpy arrays for compatibility with sklearn
    y_pred_np = y_pred.cpu().numpy().astype(int).flatten()
    y_test_np = y_test.cpu().numpy().astype(int).flatten()

    time_end = clock()  
    if SAVE_NAME != "":
        torch.save(model.state_dict(), SAVE_NAME)

if __name__ == "__main__":
    main()
