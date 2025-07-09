import torch
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
import seaborn as sns
from tqdm import tqdm  # Add tqdm import

# Set random seed for reproducibility
torch.manual_seed(42)
np.random.seed(42)

def model_train(model, train_loader, val_loader, criterion, optimizer, device, num_epochs=100, patience=10):

    model.to(device)
    
    # Initialize variables for tracking training
    best_val_loss = float('inf')
    best_model_state = None
    no_improve_epochs = 0
    
    # Lists to store metrics
    train_losses = []
    val_losses = []
    train_accs = []
    val_accs = []
    
    for epoch in range(num_epochs):
        # Training phase
        model.train()
        running_loss = 0.0
        train_topk_correct = 0
        train_total = 0
        k = 1  # Number of classes to consider for top-k accuracy
        # Add tqdm progress bar for training loop
        for inputs, labels in tqdm(train_loader, desc=f"Epoch {epoch+1}/{num_epochs} - Training", leave=False):
            inputs, labels = inputs.to(device), labels.to(device)
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            running_loss += loss.item() * inputs.size(0)
            # Compute top-k accuracy
            _, topk_preds = torch.topk(outputs, k, dim=1)
            match = topk_preds.eq(labels.view(-1, 1).expand_as(topk_preds))
            train_topk_correct += match.any(dim=1).float().sum().item()
            train_total += labels.size(0)
        epoch_train_loss = running_loss / len(train_loader.dataset)
        epoch_train_acc = train_topk_correct / train_total
        train_losses.append(epoch_train_loss)
        train_accs.append(epoch_train_acc)
        # Validation phase
        model.eval()
        val_running_loss = 0.0
        val_topk_correct = 0
        val_total = 0
        with torch.no_grad():
            for inputs, labels in tqdm(val_loader, desc=f"Epoch {epoch+1}/{num_epochs} - Validation", leave=False):
                inputs, labels = inputs.to(device), labels.to(device)
                outputs = model(inputs)
                loss = criterion(outputs, labels)
                val_running_loss += loss.item() * inputs.size(0)
                _, topk_preds = torch.topk(outputs, k, dim=1)
                match = topk_preds.eq(labels.view(-1, 1).expand_as(topk_preds))
                val_topk_correct += match.any(dim=1).float().sum().item()
                val_total += labels.size(0)
        epoch_val_loss = val_running_loss / len(val_loader.dataset)
        epoch_val_acc = val_topk_correct / val_total
        val_losses.append(epoch_val_loss)
        val_accs.append(epoch_val_acc)
        # Print progress
        print(f'Epoch {epoch+1}/{num_epochs}, '
              f'Train Loss: {epoch_train_loss:.4f}, Train Top-{k} Acc: {epoch_train_acc:.4f}, '
              f'Val Loss: {epoch_val_loss:.4f}, Val Top-{k} Acc: {epoch_val_acc:.4f}')
        # Check early stopping condition
        if epoch_val_loss < best_val_loss:
            best_val_loss = epoch_val_loss
            best_model_state = model.state_dict().copy()
            no_improve_epochs = 0
        else:
            no_improve_epochs += 1
        if no_improve_epochs >= patience:
            print(f'Early stopping at epoch {epoch+1}')
            break
    # Load best model
    model.load_state_dict(best_model_state)
    return model, {"train_losses": train_losses, "val_losses": val_losses, "train_accs": train_accs, "val_accs": val_accs}


def evaluate_model(model, test_loader, criterion, device, num_classes):
    model.eval()
    running_loss = 0.0
    all_targets = []
    all_topk_preds = []
    with torch.no_grad():
        for inputs, labels in test_loader:
            inputs, labels = inputs.to(device), labels.to(device)
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            running_loss += loss.item() * inputs.size(0)
            _, topk_preds = torch.topk(outputs, num_classes, dim=1)
            all_topk_preds.append(topk_preds.cpu().numpy())
            all_targets.extend(labels.cpu().numpy())
    # Flatten predictions
    all_topk_preds = np.concatenate(all_topk_preds, axis=0)
    for k in range(1, num_classes + 1):
        # Compute top-k accuracy: at least one of the predicted classes is correct
        correct = 0
        for i, label in enumerate(all_targets):
            if label in all_topk_preds[i][:k]:
                correct += 1
        test_acc = correct / len(all_targets)
        print(f'Test Top-{k} Accuracy: {test_acc:.4f}')
    test_loss = running_loss / len(test_loader.dataset)
    print(f'Test Loss: {test_loss:.4f}')
    # For reporting, use the top-1 prediction for confusion matrix and report
    all_preds = all_topk_preds[:, 0]
    print("\nClassification Report (Top-1):")
    print(classification_report(all_targets, all_preds, zero_division=0))
    # Get datasets object to access label mapping if needed
    dataset = test_loader.dataset
    while hasattr(dataset, 'dataset'):
        dataset = dataset.dataset
    if hasattr(dataset, 'inverse_mapping'):
        print("\nClass Index to Original Label Mapping:")
        for model_idx, original_label in dataset.inverse_mapping.items():
            print(f"  Model class {model_idx} → Original label {original_label}")
    cm = confusion_matrix(all_targets, all_preds)
    plt.figure(figsize=(12, 10))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
    plt.title('Confusion Matrix (Top-1)')
    plt.ylabel('True Label')
    plt.xlabel('Predicted Label')
    plt.savefig('confusion_matrix.png')
    return test_acc, all_preds, all_targets

