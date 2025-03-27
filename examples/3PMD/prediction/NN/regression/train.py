import torch
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
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
    
    for epoch in range(num_epochs):
        # Training phase
        model.train()
        running_loss = 0.0
        train_preds = []
        train_targets = []
        
        # Add tqdm progress bar for training loop
        for inputs, targets in tqdm(train_loader, desc=f"Epoch {epoch+1}/{num_epochs} - Training", leave=False):
            inputs, targets = inputs.to(device), targets.to(device)   
            
            optimizer.zero_grad()
            
            outputs = model(inputs)
            loss = criterion(outputs, targets)
            loss.backward()
            optimizer.step()
            
            running_loss += loss.item() * inputs.size(0)
            train_preds.extend(outputs.detach().cpu().numpy())
            train_targets.extend(targets.cpu().numpy())
            
        epoch_train_loss = running_loss / len(train_loader.dataset)
        train_losses.append(epoch_train_loss)
        
        # Validation phase
        model.eval()
        val_running_loss = 0.0
        val_preds = []
        val_targets = []
        
        # Add tqdm progress bar for validation loop
        with torch.no_grad():
            for inputs, targets in tqdm(val_loader, desc=f"Epoch {epoch+1}/{num_epochs} - Validation", leave=False):
                inputs, targets = inputs.to(device), targets.to(device)
                
                outputs = model(inputs)
                loss = criterion(outputs, targets)
                
                val_running_loss += loss.item() * inputs.size(0)
                val_preds.extend(outputs.cpu().numpy())
                val_targets.extend(targets.cpu().numpy())
        
        epoch_val_loss = val_running_loss / len(val_loader.dataset)
        val_losses.append(epoch_val_loss)
        
        # Convert predictions and targets to numpy arrays for metric calculation
        train_preds_np = np.array(train_preds)
        train_targets_np = np.array(train_targets)
        val_preds_np = np.array(val_preds)
        val_targets_np = np.array(val_targets)
        
        # Calculate MSE for individual components (sin and cos)
        train_mse_sin = mean_squared_error(train_targets_np[:, 0], train_preds_np[:, 0])
        train_mse_cos = mean_squared_error(train_targets_np[:, 1], train_preds_np[:, 1])
        val_mse_sin = mean_squared_error(val_targets_np[:, 0], val_preds_np[:, 0])
        val_mse_cos = mean_squared_error(val_targets_np[:, 1], val_preds_np[:, 1])
        
        # Print progress
        print(f'Epoch {epoch+1}/{num_epochs}, '
            f'Train Loss: {epoch_train_loss:.4f}, '
            f'Train MSE (sin): {train_mse_sin:.4f}, Train MSE (cos): {train_mse_cos:.4f}, '
            f'Val Loss: {epoch_val_loss:.4f}, '
            f'Val MSE (sin): {val_mse_sin:.4f}, Val MSE (cos): {val_mse_cos:.4f}')

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
    
    return model, {"train_losses": train_losses, "val_losses": val_losses}

def evaluate_model(model, test_loader, criterion, device):
    
    model.eval()
    running_loss = 0.0
    all_preds = []
    all_targets = []
    
    with torch.no_grad():
        for inputs, targets in test_loader:
            inputs, targets = inputs.to(device), targets.to(device)
            
            outputs = model(inputs)
            loss = criterion(outputs, targets)
            
            running_loss += loss.item() * inputs.size(0)
            all_preds.extend(outputs.cpu().numpy())
            all_targets.extend(targets.cpu().numpy())
            
    test_loss = running_loss / len(test_loader.dataset)
    
    # Convert to numpy arrays
    all_preds = np.array(all_preds)
    all_targets = np.array(all_targets)
    
    # Calculate metrics
    mse_sin = mean_squared_error(all_targets[:, 0], all_preds[:, 0])
    mse_cos = mean_squared_error(all_targets[:, 1], all_preds[:, 1])
    mae_sin = mean_absolute_error(all_targets[:, 0], all_preds[:, 0])
    mae_cos = mean_absolute_error(all_targets[:, 1], all_preds[:, 1])
    r2_sin = r2_score(all_targets[:, 0], all_preds[:, 0])
    r2_cos = r2_score(all_targets[:, 1], all_preds[:, 1])
    
    print(f'Test Loss: {test_loss:.4f}')
    print(f'MSE - Sin: {mse_sin:.4f}, Cos: {mse_cos:.4f}')
    print(f'MAE - Sin: {mae_sin:.4f}, Cos: {mae_cos:.4f}')
    print(f'R² Score - Sin: {r2_sin:.4f}, Cos: {r2_cos:.4f}')
    
    # Calculate the angular error (angle between predicted and true vectors)
    # Convert sin/cos predictions back to angles
    true_angles = np.arctan2(all_targets[:, 0], all_targets[:, 1])
    pred_angles = np.arctan2(all_preds[:, 0], all_preds[:, 1])
    
    # Calculate error in radians, accounting for periodicity
    angle_errors = np.abs(np.arctan2(
        np.sin(true_angles - pred_angles),
        np.cos(true_angles - pred_angles)
    ))
    
    mean_angle_error = np.mean(angle_errors)
    median_angle_error = np.median(angle_errors)
    
    print(f'Mean Angular Error: {mean_angle_error:.4f} radians ({np.degrees(mean_angle_error):.2f} degrees)')
    print(f'Median Angular Error: {median_angle_error:.4f} radians ({np.degrees(median_angle_error):.2f} degrees)')
    
    # Plot actual vs predicted values
    plt.figure(figsize=(12, 5))
    
    # Plot for sine
    plt.subplot(1, 2, 1)
    plt.scatter(all_targets[:, 0], all_preds[:, 0], alpha=0.5)
    plt.plot([-1, 1], [-1, 1], 'r--')
    plt.title('True vs Predicted Sin(θ)')
    plt.xlabel('True Sin(θ)')
    plt.ylabel('Predicted Sin(θ)')
    plt.grid(True)
    
    # Plot for cosine
    plt.subplot(1, 2, 2)
    plt.scatter(all_targets[:, 1], all_preds[:, 1], alpha=0.5)
    plt.plot([-1, 1], [-1, 1], 'r--')
    plt.title('True vs Predicted Cos(θ)')
    plt.xlabel('True Cos(θ)')
    plt.ylabel('Predicted Cos(θ)')
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig('regression_predictions.png')
    
    # Plot histogram of angular errors
    plt.figure(figsize=(10, 6))
    plt.hist(np.degrees(angle_errors), bins=50, alpha=0.7)
    plt.axvline(np.degrees(mean_angle_error), color='r', linestyle='--', 
                label=f'Mean Error: {np.degrees(mean_angle_error):.2f}°')
    plt.axvline(np.degrees(median_angle_error), color='g', linestyle='--', 
                label=f'Median Error: {np.degrees(median_angle_error):.2f}°')
    plt.title('Histogram of Angular Errors')
    plt.xlabel('Error (degrees)')
    plt.ylabel('Frequency')
    plt.legend()
    plt.grid(True)
    plt.savefig('angular_errors.png')
    
    return test_loss, all_preds, all_targets
