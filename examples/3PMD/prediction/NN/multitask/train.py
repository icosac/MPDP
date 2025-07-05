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
from model import MultiTaskNeuralNet  # Assuming you've renamed the model to MultiTaskNeuralNet
from tqdm import tqdm  # Progress bar

from pathlib import Path

PROJECT_PATH = os.path.dirname(os.path.abspath(__file__))
DATASETS_PATH = os.path.join(Path(PROJECT_PATH).parent, "..", "datasets")
DATASET_NAME = "big_smaller.csv"
DATASET_PATH = os.path.join(DATASETS_PATH, DATASET_NAME)
SAVE_PATH = os.path.join(PROJECT_PATH, "models")
SAVE_NAME = os.path.join(SAVE_PATH, "model.pt")
DO_PLOTS = True

# Set device to GPU if available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# device = "cpu"
print("Using device:", device)

def model_train(model, X_train, y_class_train, y_reg_train, X_val, y_class_val, y_reg_val):
    model.to(device)

    # Loss functions for both tasks
    class_loss_fn = nn.CrossEntropyLoss()  # For multi-class classification
    reg_loss_fn = nn.MSELoss()  # For regression

    optimizer = optim.AdamW(model.parameters(), lr=0.001)

    n_epochs = 100
    batch_size = 16

    train_losses, val_losses = [], []
    train_accuracies, val_accuracies = [], []

    train_loader = DataLoader(TensorDataset(X_train, y_class_train, y_reg_train), batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(TensorDataset(X_val, y_class_val, y_reg_val), batch_size=batch_size, shuffle=False)

    best_val_loss = float('inf')
    best_weights = None

    for epoch in range(n_epochs):
        model.train()
        epoch_class_loss, epoch_reg_loss = 0.0, 0.0
        correct, total = 0, 0

        progress_bar = tqdm(train_loader, desc=f"Epoch {epoch+1}/{n_epochs}", dynamic_ncols=True)

        for X_batch, y_class_batch, y_reg_batch in progress_bar:
            X_batch, y_class_batch, y_reg_batch = X_batch.to(device), y_class_batch.to(device), y_reg_batch.to(device)

            # Convert one-hot to class indices for classification
            y_class_batch = torch.argmax(y_class_batch, dim=1)

            # Reshape regression target to match the shape of the regression predictions
            y_reg_batch = y_reg_batch.view(-1, 1)  # Reshaping to [batch_size, 1]

            optimizer.zero_grad()
            # Forward pass
            class_pred, reg_pred = model(X_batch)

            # Calculate losses
            class_loss = class_loss_fn(class_pred, y_class_batch)
            reg_loss = reg_loss_fn(reg_pred, y_reg_batch)
            total_loss = class_loss + reg_loss  # Combine both losses

            total_loss.backward()
            optimizer.step()

            # Accumulate training losses and accuracy
            epoch_class_loss += class_loss.item() * X_batch.size(0)
            epoch_reg_loss += reg_loss.item() * X_batch.size(0)

            preds = torch.argmax(class_pred, dim=1)
            correct += (preds == y_class_batch).sum().item()
            total += y_class_batch.size(0)

            progress_bar.set_postfix(loss=total_loss.item())

        train_losses.append((epoch_class_loss + epoch_reg_loss) / total)
        train_accuracies.append(correct / total)

        # Validation phase
        model.eval()
        val_loss, correct, total = 0.0, 0, 0

        with torch.no_grad():
            for X_batch, y_class_batch, y_reg_batch in val_loader:
                X_batch, y_class_batch, y_reg_batch = X_batch.to(device), y_class_batch.to(device), y_reg_batch.to(device)

                y_class_batch = torch.argmax(y_class_batch, dim=1)  # Convert one-hot to class indices
                y_reg_batch = y_reg_batch.view(-1, 1)  # Reshaping to [batch_size, 1]

                class_pred, reg_pred = model(X_batch)

                class_loss = class_loss_fn(class_pred, y_class_batch)
                reg_loss = reg_loss_fn(reg_pred, y_reg_batch)
                val_loss += (class_loss + reg_loss).item() * X_batch.size(0)

                preds = torch.argmax(class_pred, dim=1)
                correct += (preds == y_class_batch).sum().item()
                total += y_class_batch.size(0)

        val_losses.append(val_loss / total)
        val_accuracies.append(correct / total)

        # Save best model based on validation loss (or accuracy)
        if val_losses[-1] < best_val_loss:
            best_val_loss = val_losses[-1]
            best_weights = copy.deepcopy(model.state_dict())

    # Restore best model
    model.load_state_dict(best_weights)

    # Plotting
    if DO_PLOTS:
        plt.figure(figsize=(12, 6))

        # Plot Loss
        plt.subplot(1, 2, 1)
        plt.plot(train_losses, label="Train Loss")
        plt.plot(val_losses, label="Validation Loss")
        plt.xlabel("Epoch")
        plt.ylabel("Loss")
        plt.legend()

        # Plot Accuracy
        plt.subplot(1, 2, 2)
        plt.plot(train_accuracies, label="Train Accuracy")
        plt.plot(val_accuracies, label="Validation Accuracy")
        plt.xlabel("Epoch")
        plt.ylabel("Accuracy")
        plt.legend()

        # plt.show()
        plt.savefig("loss_accuracy_multitask.png")

    return best_val_loss

def evaluate_model(model, X_test, y_test):
    model.eval()
    with torch.no_grad():
        X_test = X_test.to(device)
        y_test_class = torch.argmax(y_test, dim=1).to(device)  # Get class indices
        y_pred_class, y_pred_reg = model(X_test)

        y_pred_class = torch.argmax(y_pred_class, dim=1)  # Get predicted class indices
        y_pred_reg = y_pred_reg.cpu().numpy()  # Regression output

    return y_pred_class.cpu().numpy(), y_test_class.cpu().numpy(), y_pred_reg, y_test.cpu().numpy()


def main():
    assert os.path.exists(DATASET_PATH), f"Dataset {DATASET_PATH} does not exist"
    if SAVE_NAME:
        os.makedirs(SAVE_PATH, exist_ok=True)

    time_start = clock()

    data = pd.read_csv(DATASET_PATH, sep='\s+')
    features = ['kmax', 'theta_i', 'theta_f', 'alpha_m', 'alpha_f']
    target_class = 'id_man_comb'  # Classification target
    target_reg = 'th_m'  # Define your regression target column

    X = data[features].values
    y_class = data[target_class].values
    y_reg = data[target_reg].values

    # Encoding target labels for classification
    encoder = OneHotEncoder(sparse_output=False)
    y_class = encoder.fit_transform(y_class.reshape(-1, 1))

    X = torch.tensor(X, dtype=torch.float32)
    y_class = torch.tensor(y_class, dtype=torch.float32)
    y_reg = torch.tensor(y_reg, dtype=torch.float32)

    # Split data
    X_train, X_test, y_class_train, y_class_test, y_reg_train, y_reg_test = train_test_split(
        X, y_class, y_reg, test_size=0.2, stratify=y_class, random_state=42
    )

    # Train model
    model = MultiTaskNeuralNet()
    best_val_loss = model_train(model, X_train, y_class_train, y_reg_train, X_test, y_class_test, y_reg_test)

    # Evaluate model
    y_pred_class, y_test_class, y_pred_reg, y_test_reg = evaluate_model(model, X_test, y_class_test)

    # Print metrics
    print("Classification Accuracy:", accuracy_score(y_test_class, y_pred_class))
    print("Classification Precision:", precision_score(y_test_class, y_pred_class, average='weighted'))
    print("Classification Recall:", recall_score(y_test_class, y_pred_class, average='weighted'))
    print("Classification F1 Score:", f1_score(y_test_class, y_pred_class, average='weighted'))

    # You can evaluate regression metrics if needed (e.g. RMSE, R^2)
    print("Regression RMSE:", np.sqrt(np.mean((y_test_reg - y_pred_reg) ** 2)))

    # Save model
    if SAVE_NAME:
        torch.save(model.state_dict(), SAVE_NAME)


if __name__ == "__main__":
    main()
