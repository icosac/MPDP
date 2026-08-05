import os

import matplotlib.pyplot as plt
import numpy as np
import torch
from sklearn.metrics import ConfusionMatrixDisplay, accuracy_score, classification_report, confusion_matrix
from tqdm import tqdm


def model_train(model, train_loader, val_loader, criterion, optimizer, device, num_epochs=100, patience=10):
    model.to(device)

    best_val_loss = float("inf")
    best_model_state = None
    no_improve_epochs = 0
    history = {
        "train_losses": [],
        "val_losses": [],
        "train_accs": [],
        "val_accs": [],
    }

    for epoch in range(num_epochs):
        model.train()
        train_loss, train_targets, train_preds = _run_epoch(
            model, train_loader, criterion, device, optimizer=optimizer, desc=f"Epoch {epoch + 1}/{num_epochs}"
        )
        model.eval()
        val_loss, val_targets, val_preds = _run_epoch(model, val_loader, criterion, device)

        train_acc = accuracy_score(train_targets, train_preds)
        val_acc = accuracy_score(val_targets, val_preds)
        history["train_losses"].append(train_loss)
        history["val_losses"].append(val_loss)
        history["train_accs"].append(train_acc)
        history["val_accs"].append(val_acc)

        print(
            f"Epoch {epoch + 1}/{num_epochs}, "
            f"Train Loss: {train_loss:.6f}, Train Acc: {train_acc:.6f}, "
            f"Val Loss: {val_loss:.6f}, Val Acc: {val_acc:.6f}"
        )

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_model_state = {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}
            no_improve_epochs = 0
        else:
            no_improve_epochs += 1

        if no_improve_epochs >= patience:
            print(f"Early stopping at epoch {epoch + 1}")
            break

    if best_model_state is not None:
        model.load_state_dict(best_model_state)

    return model, history


def evaluate_model(model, data_loader, criterion, device, num_classes, plot_path=None, split_name="validation"):
    model.eval()
    loss, targets, preds, topk_preds = _predict(model, data_loader, criterion, device, num_classes)

    print(f"{split_name.capitalize()} Loss: {loss:.6f}")
    for k in range(1, num_classes + 1):
        correct = sum(target in topk_preds[idx, :k] for idx, target in enumerate(targets))
        print(f"{split_name.capitalize()} Top-{k} Accuracy: {correct / len(targets):.6f}")

    print("\nClassification Report:")
    print(classification_report(targets, preds, zero_division=0))

    dataset = data_loader.dataset
    if hasattr(dataset, "inverse_mapping"):
        print("\nClass Index to Original Label Mapping:")
        for model_idx, original_label in sorted(dataset.inverse_mapping.items()):
            print(f"  Model class {model_idx}: original man {original_label}")

    if plot_path:
        os.makedirs(plot_path, exist_ok=True)
        cm = confusion_matrix(targets, preds, normalize="true")
        fig, ax = plt.subplots(figsize=(12, 10))
        disp = ConfusionMatrixDisplay(confusion_matrix=cm)
        disp.plot(cmap="Blues", values_format=".4f", xticks_rotation=45, ax=ax)
        fig.tight_layout()
        fig.savefig(os.path.join(plot_path, f"confusion_matrix_norm_{split_name}.png"), dpi=300)
        plt.close(fig)
        np.savetxt(
            os.path.join(plot_path, f"confusion_matrix_norm_{split_name}.csv"),
            cm,
            delimiter=",",
            fmt="%.6f",
        )

    return accuracy_score(targets, preds), preds, targets


def plot_training_results(history, plot_path):
    os.makedirs(plot_path, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].plot(history["train_losses"], label="Train Loss")
    axes[0].plot(history["val_losses"], label="Validation Loss")
    axes[0].set_title("Loss over Epochs")
    axes[0].set_xlabel("Epoch")
    axes[0].set_ylabel("Loss")
    axes[0].legend()

    axes[1].plot(history["train_accs"], label="Train Accuracy")
    axes[1].plot(history["val_accs"], label="Validation Accuracy")
    axes[1].set_title("Accuracy over Epochs")
    axes[1].set_xlabel("Epoch")
    axes[1].set_ylabel("Accuracy")
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(os.path.join(plot_path, "training_history.png"), dpi=300)
    plt.close(fig)


def _run_epoch(model, data_loader, criterion, device, optimizer=None, desc=None):
    running_loss = 0.0
    targets = []
    preds = []
    iterable = tqdm(data_loader, desc=desc, leave=False) if desc else data_loader

    with torch.set_grad_enabled(optimizer is not None):
        for inputs, labels in iterable:
            inputs = inputs.to(device)
            labels = labels.to(device)

            if optimizer is not None:
                optimizer.zero_grad()

            outputs = model(inputs)
            loss = criterion(outputs, labels)

            if optimizer is not None:
                loss.backward()
                optimizer.step()

            running_loss += loss.item() * inputs.size(0)
            preds.extend(torch.argmax(outputs, dim=1).cpu().numpy())
            targets.extend(labels.cpu().numpy())

    return running_loss / len(data_loader.dataset), np.array(targets), np.array(preds)


def _predict(model, data_loader, criterion, device, num_classes):
    running_loss = 0.0
    targets = []
    preds = []
    topk = []

    with torch.no_grad():
        for inputs, labels in data_loader:
            inputs = inputs.to(device)
            labels = labels.to(device)
            outputs = model(inputs)
            loss = criterion(outputs, labels)

            running_loss += loss.item() * inputs.size(0)
            preds.extend(torch.argmax(outputs, dim=1).cpu().numpy())
            targets.extend(labels.cpu().numpy())
            topk.append(torch.topk(outputs, num_classes, dim=1).indices.cpu().numpy())

    return running_loss / len(data_loader.dataset), np.array(targets), np.array(preds), np.concatenate(topk, axis=0)
