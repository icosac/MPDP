import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, random_split
import matplotlib.pyplot as plt
import numpy as np
from model import NeuralNet
from dataset import DubinsDataset
from train import model_train, evaluate_model

DO_TRAINING = True

#################################################
################ PLOT FUNCS #####################
#################################################

def plot_training_results(history):
    """
    Plots training and validation metrics.
    
    Args:
        history: Dictionary containing training history
    """
    plt.figure(figsize=(10, 5))
    
    # Plot loss
    plt.plot(history['train_losses'], label='Train Loss')
    plt.plot(history['val_losses'], label='Validation Loss')
    plt.title('Loss over Epochs')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('training_history.png')

#################################################
################ INFERENCE CLASS ################
#################################################

class ModelInference:
    def __init__(self, model_path, scaler, input_size, hidden_size, device='cpu'):
        """
        Class for inference with a trained model.
        
        Args:
            model_path: Path to the saved model
            scaler: StandardScaler fit on training data
            input_size: Number of input features
            hidden_size: Number of hidden units
            device: Device to use for computation
        """
        self.device = device
        self.scaler = scaler
        
        # Initialize model (output size is 2 for sin and cos)
        self.model = NeuralNet(input_size, hidden_size, out_size=2)
        self.model.load_state_dict(torch.load(model_path, map_location=device))
        self.model.to(device)
        self.model.eval()
        
    def predict(self, features):
        """
        Predicts sin and cos values for the angle.
        
        Args:
            features: Features as numpy array shape (5,) or (n, 5)
        
        Returns:
            Tuple of (sin_th_m, cos_th_m, predicted_angle_radians, predicted_angle_degrees)
        """
        # Handle single sample vs batch
        single_sample = False
        if len(features.shape) == 1:
            features = features.reshape(1, -1)
            single_sample = True
            
        # Preprocess features
        features = self.scaler.transform(features)
        features = torch.FloatTensor(features).to(self.device)
        
        # Get prediction
        with torch.no_grad():
            outputs = self.model(features)
            
        # Convert to numpy
        sin_cos_values = outputs.cpu().numpy()
        
        # Calculate angles
        predicted_angles_rad = np.arctan2(sin_cos_values[:, 0], sin_cos_values[:, 1])
        predicted_angles_deg = np.degrees(predicted_angles_rad)
        
        if single_sample:
            return (sin_cos_values[0, 0], sin_cos_values[0, 1], 
                    predicted_angles_rad[0], predicted_angles_deg[0])
        else:
            return (sin_cos_values[:, 0], sin_cos_values[:, 1], 
                    predicted_angles_rad, predicted_angles_deg)

#################################################
################ MAIN FUNCTION ##################
#################################################

def main(data_path, model_save_path='model_regression.pt'):
    
    try:
        with open(data_path, 'r') as f:
            print("Data preview:")
            for i, line in enumerate(f):
                print(line.strip())
                if i >= 5:  # Show first 5 lines
                    break
    except Exception as e:
        print(f"Error reading file: {e}")
    
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # Load dataset
    dataset = DubinsDataset(data_path)
    print(f"Dataset loaded: {len(dataset)} samples")
    
    # Split dataset
    train_size = int(0.7 * len(dataset))
    val_size = int(0.15 * len(dataset))
    test_size = len(dataset) - train_size - val_size
    
    train_dataset, val_dataset, test_dataset = random_split(
        dataset, [train_size, val_size, test_size]
    )
    
    # Create data loaders
    batch_size = 32
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size)
    test_loader = DataLoader(test_dataset, batch_size=batch_size)
    
    # Initialize model
    input_size = 5  # Number of features
    hidden_size = 64
    output_size = 2  # sin and cos components
    
    model = NeuralNet(input_size, hidden_size, output_size)
    
    # Define loss function and optimizer
    criterion = nn.MSELoss()  # Mean Squared Error loss for regression
    optimizer = optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-5)
    
    if DO_TRAINING:
        # Train model       
        print("Starting training...")
        trained_model, history = model_train(
            model, train_loader, val_loader, criterion, optimizer, device, 
            num_epochs=100, patience=10
        )
        
        # Save the model
        torch.save(trained_model.state_dict(), model_save_path)
        print(f"Model saved to {model_save_path}")
        
        # Plot training results
        plot_training_results(history)
        
        # Evaluate on test set
        print("\nEvaluating on test set...")
        test_loss, _, _ = evaluate_model(trained_model, test_loader, criterion, device)
    
    scaler = dataset.get_scaler()
    inference = ModelInference(model_save_path, scaler, input_size, hidden_size, device)
    
    # Example features for prediction
    example_features = np.array([1, 2.3562, 2.3562, 1.5708, 2.3562])
    
    # Predict
    sin_val, cos_val, angle_rad, angle_deg = inference.predict(example_features)
    print(f"Example features: {example_features}")
    print(f"Predicted sin(th_m): {sin_val:.4f}, cos(th_m): {cos_val:.4f}")
    print(f"Predicted angle: {angle_rad:.4f} radians ({angle_deg:.2f} degrees)")
    
    
if __name__ == "__main__":
    main('/home/davide/Desktop/MPDP/examples/3PMD/prediction/NN/regression/datasets/small.csv')