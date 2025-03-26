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
    # Create a figure with 2 rows and 2 columns
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot combined loss
    axs[0, 0].plot(history['train_losses'], label='Train Loss')
    axs[0, 0].plot(history['val_losses'], label='Validation Loss')
    axs[0, 0].set_title('Combined Loss over Epochs')
    axs[0, 0].set_xlabel('Epoch')
    axs[0, 0].set_ylabel('Loss')
    axs[0, 0].legend()
    axs[0, 0].grid(True)
    
    # Plot classification loss
    axs[0, 1].plot(history['train_class_losses'], label='Train Classification Loss')
    axs[0, 1].plot(history['val_class_losses'], label='Validation Classification Loss')
    axs[0, 1].set_title('Classification Loss over Epochs')
    axs[0, 1].set_xlabel('Epoch')
    axs[0, 1].set_ylabel('Loss')
    axs[0, 1].legend()
    axs[0, 1].grid(True)
    
    # Plot regression loss
    axs[1, 0].plot(history['train_reg_losses'], label='Train Regression Loss')
    axs[1, 0].plot(history['val_reg_losses'], label='Validation Regression Loss')
    axs[1, 0].set_title('Regression Loss over Epochs')
    axs[1, 0].set_xlabel('Epoch')
    axs[1, 0].set_ylabel('Loss')
    axs[1, 0].legend()
    axs[1, 0].grid(True)
    
    # Remove the unused subplot
    fig.delaxes(axs[1, 1])
    
    plt.tight_layout()
    plt.savefig('training_history.png')
    plt.close()

#################################################
################ INFERENCE CLASS ################
#################################################

class ModelInference:
    def __init__(self, model_path, scaler, input_size, hidden_size, num_classes, device='cpu'):
        """
        Class for inference with a trained multi-task model.
        
        Args:
            model_path: Path to the saved model
            scaler: StandardScaler fit on training data
            input_size: Number of input features
            hidden_size: Number of hidden units
            num_classes: Number of classification classes
            device: Device to use for computation
        """
        self.device = device
        self.scaler = scaler
        
        # Initialize model
        self.model = NeuralNet(input_size, hidden_size, class_out_size=num_classes, reg_out_size=2)
        self.model.load_state_dict(torch.load(model_path, map_location=device))
        self.model.to(device)
        self.model.eval()
        
        # Store the label mapping
        self.label_mapping = None
        self.inverse_mapping = None
        
    def set_label_mapping(self, label_mapping, inverse_mapping):
        """
        Set the label mapping for classification results.
        
        Args:
            label_mapping: Dictionary mapping original labels to indices
            inverse_mapping: Dictionary mapping indices to original labels
        """
        self.label_mapping = label_mapping
        self.inverse_mapping = inverse_mapping
        
    def predict(self, features):
        """
        Predicts class and angle for input features.
        
        Args:
            features: Features as numpy array shape (5,) or (n, 5)
        
        Returns:
            Dictionary containing class prediction and angle prediction
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
            class_outputs, reg_outputs = self.model(features)
            
        # Get class predictions
        _, predicted_classes = torch.max(class_outputs, 1)
        predicted_classes = predicted_classes.cpu().numpy()
        
        # Get original class labels if mapping is available
        if self.inverse_mapping:
            original_classes = np.array([self.inverse_mapping[idx] for idx in predicted_classes])
        else:
            original_classes = predicted_classes
            
        # Get regression predictions
        sin_cos_values = reg_outputs.cpu().numpy()
        
        # Calculate angles
        predicted_angles_rad = np.arctan2(sin_cos_values[:, 0], sin_cos_values[:, 1])
        predicted_angles_deg = np.degrees(predicted_angles_rad)
        
        if single_sample:
            return {
                'class_idx': predicted_classes[0],
                'class_original': original_classes[0],
                'sin_th_m': sin_cos_values[0, 0],
                'cos_th_m': sin_cos_values[0, 1],
                'angle_rad': predicted_angles_rad[0],
                'angle_deg': predicted_angles_deg[0]
            }
        else:
            return {
                'class_idx': predicted_classes,
                'class_original': original_classes,
                'sin_th_m': sin_cos_values[:, 0],
                'cos_th_m': sin_cos_values[:, 1],
                'angle_rad': predicted_angles_rad,
                'angle_deg': predicted_angles_deg
            }

#################################################
################ MAIN FUNCTION ##################
#################################################

def main(data_path, model_save_path='model_multitask.pt'):
    
    # Preview data
    try:
        with open(data_path, 'r') as f:
            print("Data preview:")
            for i, line in enumerate(f):
                print(line.strip())
                if i >= 5:  # Show first 5 lines
                    break
    except Exception as e:
        print(f"Error reading file: {e}")
    
    # Set device
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # Load dataset
    dataset = DubinsDataset(data_path)
    print(f"Dataset loaded: {len(dataset)} samples")
    print(f"Number of classes: {dataset.num_classes}")
    
    # Get label mapping for evaluation
    label_mapping, inverse_mapping = dataset.get_label_mapping()
    
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
    num_classes = dataset.num_classes  # Number of classes for classification
    reg_output_size = 2  # sin and cos components
    
    model = NeuralNet(input_size, hidden_size, class_out_size=num_classes, reg_out_size=reg_output_size)
    
    # Define loss functions and optimizer
    criterion_class = nn.CrossEntropyLoss()  # Cross entropy loss for classification
    criterion_reg = nn.MSELoss()  # Mean squared error loss for regression
    optimizer = optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-5)
    
    if DO_TRAINING:
        # Train model       
        print("Starting training...")
        trained_model, history = model_train(
            model, train_loader, val_loader, criterion_class, criterion_reg, optimizer, device, 
            num_epochs=100, patience=10
        )
        
        # Save the model
        torch.save(trained_model.state_dict(), model_save_path)
        print(f"Model saved to {model_save_path}")
        
        # Plot training results
        plot_training_results(history)
        
        # Evaluate on test set
        print("\nEvaluating on test set...")
        test_loss, metrics = evaluate_model(
            trained_model, test_loader, criterion_class, criterion_reg, device, inverse_mapping
        )
    else:
        # Load trained model
        model.load_state_dict(torch.load(model_save_path, map_location=device))
        
        # Evaluate on test set
        print("\nEvaluating pre-trained model on test set...")
        test_loss, metrics = evaluate_model(
            model, test_loader, criterion_class, criterion_reg, device, inverse_mapping
        )
    
    # Initialize inference class
    scaler = dataset.get_scaler()
    inference = ModelInference(model_save_path, scaler, input_size, hidden_size, num_classes, device)
    inference.set_label_mapping(label_mapping, inverse_mapping)
    
    # Example features for prediction
    example_features = np.array([1, 2.3562, 2.3562, 1.5708, 2.3562])
    
    # Predict
    prediction = inference.predict(example_features)
    print("\nExample prediction:")
    print(f"Input features: {example_features}")
    print(f"Predicted class index: {prediction['class_idx']}")
    print(f"Predicted original class: {prediction['class_original']}")
    print(f"Predicted sin(th_m): {prediction['sin_th_m']:.4f}, cos(th_m): {prediction['cos_th_m']:.4f}")
    print(f"Predicted angle: {prediction['angle_rad']:.4f} radians ({prediction['angle_deg']:.2f} degrees)")
    
    
if __name__ == "__main__":
    main('/home/davide/Desktop/MPDP/examples/3PMD/prediction/NN/regression/datasets/small.csv')