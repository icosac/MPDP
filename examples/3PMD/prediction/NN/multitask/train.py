import torch
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score, accuracy_score, confusion_matrix, classification_report
import seaborn as sns

# Set random seed for reproducibility
torch.manual_seed(42)
np.random.seed(42)

def model_train(model, train_loader, val_loader, criterion_class, criterion_reg, optimizer, device, num_epochs=100, patience=10):
    """
    Train a multi-task model for both classification and regression.
    
    Args:
        model: The neural network model
        train_loader: DataLoader for training data
        val_loader: DataLoader for validation data
        criterion_class: Loss function for classification
        criterion_reg: Loss function for regression
        optimizer: Optimizer for training
        device: Device to use for training
        num_epochs: Maximum number of epochs
        patience: Number of epochs to wait before early stopping
        
    Returns:
        Trained model and training history
    """
    model.to(device)
    
    # Initialize variables for tracking training
    best_val_loss = float('inf')
    best_model_state = None
    no_improve_epochs = 0
    
    # Lists to store metrics
    train_losses = []
    val_losses = []
    train_class_losses = []
    train_reg_losses = []
    val_class_losses = []
    val_reg_losses = []
    
    for epoch in range(num_epochs):
        # Training phase
        model.train()
        running_loss = 0.0
        running_class_loss = 0.0
        running_reg_loss = 0.0
        train_true_labels = []
        train_pred_labels = []
        train_true_targets = []
        train_pred_targets = []
            
        for inputs, targets, labels in train_loader:
            inputs, targets, labels = inputs.to(device), targets.to(device), labels.to(device)
            
            optimizer.zero_grad()
            
            # Forward pass
            out_class, out_reg = model(inputs)
            
            # Calculate losses
            loss_class = criterion_class(out_class, labels)
            loss_reg = criterion_reg(out_reg, targets)
            
            # Combined loss
            loss = loss_class + loss_reg
            
            # Backward pass
            loss.backward()
            optimizer.step()
            
            # Update running losses
            running_loss += loss.item() * inputs.size(0)
            running_class_loss += loss_class.item() * inputs.size(0)
            running_reg_loss += loss_reg.item() * inputs.size(0)
            
            # Store predictions and true values
            _, predicted_labels = torch.max(out_class, 1)
            train_true_labels.extend(labels.cpu().numpy())
            train_pred_labels.extend(predicted_labels.cpu().numpy())
            train_true_targets.extend(targets.cpu().numpy())
            train_pred_targets.extend(out_reg.detach().cpu().numpy())
            
        # Calculate epoch losses
        dataset_size = len(train_loader.dataset)
        epoch_train_loss = running_loss / dataset_size
        epoch_train_class_loss = running_class_loss / dataset_size
        epoch_train_reg_loss = running_reg_loss / dataset_size
        
        # Store losses
        train_losses.append(epoch_train_loss)
        train_class_losses.append(epoch_train_class_loss)
        train_reg_losses.append(epoch_train_reg_loss)
        
        # Calculate training metrics
        train_accuracy = accuracy_score(train_true_labels, train_pred_labels)
        train_true_targets = np.array(train_true_targets)
        train_pred_targets = np.array(train_pred_targets)
        train_mse_sin = mean_squared_error(train_true_targets[:, 0], train_pred_targets[:, 0])
        train_mse_cos = mean_squared_error(train_true_targets[:, 1], train_pred_targets[:, 1])
        
        # Validation phase
        model.eval()
        val_running_loss = 0.0
        val_running_class_loss = 0.0
        val_running_reg_loss = 0.0
        val_true_labels = []
        val_pred_labels = []
        val_true_targets = []
        val_pred_targets = []
        
        with torch.no_grad():
            for inputs, targets, labels in val_loader:
                inputs, targets, labels = inputs.to(device), targets.to(device), labels.to(device)
                
                # Forward pass
                out_class, out_reg = model(inputs)
                
                # Calculate losses
                loss_class = criterion_class(out_class, labels)
                loss_reg = criterion_reg(out_reg, targets)
                loss = loss_class + loss_reg
                
                # Update running losses
                val_running_loss += loss.item() * inputs.size(0)
                val_running_class_loss += loss_class.item() * inputs.size(0)
                val_running_reg_loss += loss_reg.item() * inputs.size(0)
                
                # Store predictions and true values
                _, predicted_labels = torch.max(out_class, 1)
                val_true_labels.extend(labels.cpu().numpy())
                val_pred_labels.extend(predicted_labels.cpu().numpy())
                val_true_targets.extend(targets.cpu().numpy())
                val_pred_targets.extend(out_reg.cpu().numpy())
        
        # Calculate epoch validation losses
        epoch_val_loss = val_running_loss / len(val_loader.dataset)
        epoch_val_class_loss = val_running_class_loss / len(val_loader.dataset)
        epoch_val_reg_loss = val_running_reg_loss / len(val_loader.dataset)
        
        # Store validation losses
        val_losses.append(epoch_val_loss)
        val_class_losses.append(epoch_val_class_loss)
        val_reg_losses.append(epoch_val_reg_loss)
        
        # Calculate validation metrics
        val_accuracy = accuracy_score(val_true_labels, val_pred_labels)
        val_true_targets = np.array(val_true_targets)
        val_pred_targets = np.array(val_pred_targets)
        val_mse_sin = mean_squared_error(val_true_targets[:, 0], val_pred_targets[:, 0])
        val_mse_cos = mean_squared_error(val_true_targets[:, 1], val_pred_targets[:, 1])
        
        # Print progress
        print(f'Epoch {epoch+1}/{num_epochs}, '
              f'Train Loss: {epoch_train_loss:.4f} (Class: {epoch_train_class_loss:.4f}, Reg: {epoch_train_reg_loss:.4f}), '
              f'Train Acc: {train_accuracy:.4f}, '
              f'Train MSE (sin/cos): {train_mse_sin:.4f}/{train_mse_cos:.4f}, '
              f'Val Loss: {epoch_val_loss:.4f} (Class: {epoch_val_class_loss:.4f}, Reg: {epoch_val_reg_loss:.4f}), '
              f'Val Acc: {val_accuracy:.4f}, '
              f'Val MSE (sin/cos): {val_mse_sin:.4f}/{val_mse_cos:.4f}')

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
    
    return model, {
        "train_losses": train_losses, 
        "val_losses": val_losses,
        "train_class_losses": train_class_losses,
        "train_reg_losses": train_reg_losses,
        "val_class_losses": val_class_losses,
        "val_reg_losses": val_reg_losses
    }

def evaluate_model(model, test_loader, criterion_class, criterion_reg, device, label_inverse_mapping=None):
    """
    Evaluate a trained multi-task model on test data.
    
    Args:
        model: The trained neural network model
        test_loader: DataLoader for test data
        criterion_class: Loss function for classification
        criterion_reg: Loss function for regression
        device: Device to use for evaluation
        label_inverse_mapping: Mapping from class indices to original labels
        
    Returns:
        test_loss, classification metrics, regression metrics
    """
    model.eval()
    running_loss = 0.0
    running_class_loss = 0.0
    running_reg_loss = 0.0
    
    all_true_labels = []
    all_pred_labels = []
    all_true_targets = []
    all_pred_targets = []
    
    with torch.no_grad():
        for inputs, targets, labels in test_loader:
            inputs, targets, labels = inputs.to(device), targets.to(device), labels.to(device)
            
            # Forward pass
            out_class, out_reg = model(inputs)
            
            # Calculate losses
            loss_class = criterion_class(out_class, labels)
            loss_reg = criterion_reg(out_reg, targets)
            loss = loss_class + loss_reg
            
            # Update running losses
            running_loss += loss.item() * inputs.size(0)
            running_class_loss += loss_class.item() * inputs.size(0)
            running_reg_loss += loss_reg.item() * inputs.size(0)
            
            # Store predictions and true values
            _, predicted_labels = torch.max(out_class, 1)
            all_true_labels.extend(labels.cpu().numpy())
            all_pred_labels.extend(predicted_labels.cpu().numpy())
            all_true_targets.extend(targets.cpu().numpy())
            all_pred_targets.extend(out_reg.cpu().numpy())
            
    # Calculate test loss
    test_loss = running_loss / len(test_loader.dataset)
    test_class_loss = running_class_loss / len(test_loader.dataset)
    test_reg_loss = running_reg_loss / len(test_loader.dataset)
    
    # Convert to numpy arrays
    all_true_labels = np.array(all_true_labels)
    all_pred_labels = np.array(all_pred_labels)
    all_true_targets = np.array(all_true_targets)
    all_pred_targets = np.array(all_pred_targets)
    
    # Calculate classification metrics
    accuracy = accuracy_score(all_true_labels, all_pred_labels)
    
    # Map back to original labels for better readability if mapping is provided
    if label_inverse_mapping:
        all_true_labels_orig = np.array([label_inverse_mapping[label] for label in all_true_labels])
        all_pred_labels_orig = np.array([label_inverse_mapping[label] for label in all_pred_labels])
        report = classification_report(all_true_labels_orig, all_pred_labels_orig)
    else:
        report = classification_report(all_true_labels, all_pred_labels)
    
    conf_matrix = confusion_matrix(all_true_labels, all_pred_labels)
    
    # Calculate regression metrics
    mse_sin = mean_squared_error(all_true_targets[:, 0], all_pred_targets[:, 0])
    mse_cos = mean_squared_error(all_true_targets[:, 1], all_pred_targets[:, 1])
    mae_sin = mean_absolute_error(all_true_targets[:, 0], all_pred_targets[:, 0])
    mae_cos = mean_absolute_error(all_true_targets[:, 1], all_pred_targets[:, 1])
    r2_sin = r2_score(all_true_targets[:, 0], all_pred_targets[:, 0])
    r2_cos = r2_score(all_true_targets[:, 1], all_pred_targets[:, 1])
    
    # Print results
    print(f'Test Loss: {test_loss:.4f} (Class: {test_class_loss:.4f}, Reg: {test_reg_loss:.4f})')
    print(f'Classification Accuracy: {accuracy:.4f}')
    print(f'Classification Report:\n{report}')
    print(f'Regression Metrics:')
    print(f'MSE - Sin: {mse_sin:.4f}, Cos: {mse_cos:.4f}')
    print(f'MAE - Sin: {mae_sin:.4f}, Cos: {mae_cos:.4f}')
    print(f'R² Score - Sin: {r2_sin:.4f}, Cos: {r2_cos:.4f}')
    
    # Calculate angular error
    true_angles = np.arctan2(all_true_targets[:, 0], all_true_targets[:, 1])
    pred_angles = np.arctan2(all_pred_targets[:, 0], all_pred_targets[:, 1])
    
    # Calculate error in radians, accounting for periodicity
    angle_errors = np.abs(np.arctan2(
        np.sin(true_angles - pred_angles),
        np.cos(true_angles - pred_angles)
    ))
    
    mean_angle_error = np.mean(angle_errors)
    median_angle_error = np.median(angle_errors)
    
    print(f'Mean Angular Error: {mean_angle_error:.4f} radians ({np.degrees(mean_angle_error):.2f} degrees)')
    print(f'Median Angular Error: {median_angle_error:.4f} radians ({np.degrees(median_angle_error):.2f} degrees)')
    
   
    return test_loss, {
        'accuracy': accuracy,
        'mse_sin': mse_sin,
        'mse_cos': mse_cos,
        'mae_sin': mae_sin,
        'mae_cos': mae_cos,
        'r2_sin': r2_sin,
        'r2_cos': r2_cos,
        'mean_angle_error': mean_angle_error,
        'median_angle_error': median_angle_error
    }