import torch
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
import seaborn as sns

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
        train_preds = []
        train_targets = []
            
        for inputs, labels in train_loader:
            inputs, labels = inputs.to(device), labels.to(device)   
            
            optimizer.zero_grad()
            
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            
            running_loss += loss.item() * inputs.size(0)
            _, pred = torch.max(outputs, 1)
            train_preds.extend(pred.cpu().numpy())
            train_targets.extend(labels.cpu().numpy())
            
        epoch_train_loss = running_loss / len(train_loader.dataset)
        epoch_train_acc = accuracy_score(train_targets, train_preds)
        train_losses.append(epoch_train_loss)
        train_accs.append(epoch_train_acc)
        
        # Validation phase
        model.eval()
        val_running_loss = 0.0
        val_preds = []
        val_targets = []
        
        with torch.no_grad():
            for inputs, labels in val_loader:
                inputs, labels = inputs.to(device), labels.to(device)
                
                outputs = model(inputs)
                loss = criterion(outputs, labels)
                
                val_running_loss += loss.item() * inputs.size(0)
                _, pred = torch.max(outputs, 1)
                val_preds.extend(pred.cpu().numpy())
                val_targets.extend(labels.cpu().numpy())
        
        epoch_val_loss = val_running_loss / len(val_loader.dataset)
        epoch_val_acc = accuracy_score(val_targets, val_preds)
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
            _, pred = torch.max(outputs, 1)
            all_preds.extend(pred.cpu().numpy())
            all_targets.extend(labels.cpu().numpy())
            
    test_loss = running_loss / len(test_loader.dataset)
    test_acc = accuracy_score(all_targets, all_preds)
    
    print(f'Test Loss: {test_loss:.4f}, Test Accuracy: {test_acc:.4f}')
    
    # Generate classification report
    print("\nClassification Report:")
    print(classification_report(all_targets, all_preds))
    
    # Get datasets object to access label mapping if needed
    dataset = test_loader.dataset
    # If this is a Subset (from random_split), get the original dataset
    while hasattr(dataset, 'dataset'):
        dataset = dataset.dataset
    
    # Print mapping from model indices to original class labels if available
    if hasattr(dataset, 'inverse_mapping'):
        print("\nClass Index to Original Label Mapping:")
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

