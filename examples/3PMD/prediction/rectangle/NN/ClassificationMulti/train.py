import torch
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
import seaborn as sns
from tqdm import tqdm  # Add tqdm import

# Set random seed for reproducibility
torch.manual_seed(42)
np.random.seed(42)

def top_k_multi_label_accuracy(y_true, y_pred_logits, k=4):
    # y_true: (batch, num_classes) multi-hot, where 1 means the class is a valid label for this sample
    # y_pred_logits: (batch, num_classes) raw logits
    # For each sample, prediction is correct if any of the top-k predicted classes is in the set of valid classes
    y_pred_topk = torch.topk(y_pred_logits, k=k, dim=1).indices.cpu().numpy()
    y_true_indices = [set(np.where(row > 0.5)[0]) for row in y_true.cpu().numpy()]
    correct = 0
    for idx, topk in enumerate(y_pred_topk):
        if any(i in y_true_indices[idx] for i in topk):
            correct += 1
    return correct / len(y_true)

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
        train_preds = []
        train_targets = []
            
        # Add tqdm progress bar for training loop
        for inputs, labels in tqdm(train_loader, desc=f"Epoch {epoch+1}/{num_epochs} - Training", leave=False):
            inputs, labels = inputs.to(device), labels.to(device)   
            
            optimizer.zero_grad()
            
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            
            running_loss += loss.item() * inputs.size(0)
            train_preds.append(outputs.detach())
            train_targets.append(labels.detach())
            
        epoch_train_loss = running_loss / len(train_loader.dataset)
        train_preds_tensor = torch.cat(train_preds, dim=0)
        train_targets_tensor = torch.cat(train_targets, dim=0)
        epoch_train_acc = top_k_multi_label_accuracy(train_targets_tensor, train_preds_tensor, k=4)
        train_losses.append(epoch_train_loss)
        train_accs.append(epoch_train_acc)
        
        # Validation phase
        model.eval()
        val_running_loss = 0.0
        val_preds = []
        val_targets = []
        
        # Add tqdm progress bar for validation loop
        with torch.no_grad():
            for inputs, labels in tqdm(val_loader, desc=f"Epoch {epoch+1}/{num_epochs} - Validation", leave=False):
                inputs, labels = inputs.to(device), labels.to(device)
                
                outputs = model(inputs)
                loss = criterion(outputs, labels)
                
                val_running_loss += loss.item() * inputs.size(0)
                val_preds.append(outputs)
                val_targets.append(labels)
        
        epoch_val_loss = val_running_loss / len(val_loader.dataset)
        val_preds_tensor = torch.cat(val_preds, dim=0)
        val_targets_tensor = torch.cat(val_targets, dim=0)
        epoch_val_acc = top_k_multi_label_accuracy(val_targets_tensor, val_preds_tensor, k=4)
        val_losses.append(epoch_val_loss)
        val_accs.append(epoch_val_acc)
        
        # Print progress
        print(f'Epoch {epoch+1}/{num_epochs}, '
            f'Train Loss: {epoch_train_loss:.4f}, Train Acc: {epoch_train_acc:.4f}, '
            f'Val Loss: {epoch_val_loss:.4f}, Val Acc: {epoch_val_acc:.4f}')

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
    
    return model, {"train_losses": train_losses, "val_losses": val_losses, 
                "train_accs": train_accs, "val_accs": val_accs}

def evaluate_model(model, test_loader, criterion, device):
    
    model.eval()
    running_loss = 0.0
    all_preds = []
    all_targets = []
    
    with torch.no_grad():
        for inputs, labels in test_loader:
            inputs, labels = inputs.to(device), labels.to(device)
            
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            
            running_loss += loss.item() * inputs.size(0)
            all_preds.append(outputs)
            all_targets.append(labels)
    test_loss = running_loss / len(test_loader.dataset)
    all_preds_tensor = torch.cat(all_preds, dim=0)
    all_targets_tensor = torch.cat(all_targets, dim=0)
    test_acc = top_k_multi_label_accuracy(all_targets_tensor, all_preds_tensor, k=4)
    print(f'Test Loss: {test_loss:.4f}, Top-4 Accuracy: {test_acc:.4f}')
    
    return test_acc, all_preds_tensor, all_targets_tensor
        for model_idx, original_label in dataset.inverse_mapping.items():
            print(f"  Model class {model_idx} → Original label {original_label}")
    
    # Plot confusion matrix
    cm = confusion_matrix(all_targets, all_preds)
    plt.figure(figsize=(12, 10))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
    plt.title('Confusion Matrix')
    plt.ylabel('True Label')
    plt.xlabel('Predicted Label')
    plt.savefig('confusion_matrix.png')
    
    return test_acc, all_preds, all_targets

