#################################################
############### REGRESSION ######################
#################################################

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, random_split
import matplotlib.pyplot as plt
import numpy as np
from model import NeuralNet
from dataset import DubinsDatasetRectangle
from train import model_train, evaluate_model

import os
from pathlib import Path

BATCH_SIZE  = 16
EPOCHS      = 2
PATIENCE    = 10
LEARN_RATE  = 0.0001
WEIGHT_DEC  = 1e-3
HIDDEN_SIZE = 128

TRIG_FUNCS   = False  # Use trigonometric features

DO_TRAINING    = True
USE_ONLY_CPU   = False
EXPORT_TO_ONNX = False

THIS_FILE_PATH = os.path.abspath(__file__)
PROJECT_PATH   = Path(THIS_FILE_PATH).parent
DATASET_PATH   = os.path.join(PROJECT_PATH.parent.parent, "datasets")
MODELS_PATH    = os.path.join(PROJECT_PATH, "models")
PLOT_PATH      = os.path.join(PROJECT_PATH, "plots")

DATASET_NAME   = os.path.join("/Users/enrico/Projects/mpdp/small_rec.csv")
MODEL_NAME     = os.path.join(MODELS_PATH, 'model_regression.pt')
ONNX_NAME      = os.path.join(MODELS_PATH, 'model_regression.onnx')

os.makedirs(MODELS_PATH, exist_ok=True)
os.makedirs(PLOT_PATH,   exist_ok=True)

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
    plt.savefig(os.path.join(PLOT_PATH, 'training_history.png'))

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
################  ONNX EXPORT  ##################
#################################################

class ScaledNeuralNet(nn.Module):
    def __init__(self, model, scaler, device='cpu'):
        super(ScaledNeuralNet, self).__init__()
        self.model = model
        self.register_buffer('mean', torch.tensor(scaler.mean_, dtype=torch.float32).to(device))
        self.register_buffer('scale', torch.tensor(scaler.scale_, dtype=torch.float32).to(device))
        
    def forward(self, x):
        x = (x - self.mean) / self.scale
        return self.model(x)

def export_to_onnx(model, save_path, input_size, scaler):
    """
    Export model to ONNX format, handling device issues correctly
    
    Args:
        model: PyTorch model to export
        save_path: Path where the ONNX model will be saved
        input_size: Size of the input tensor (number of features)
    """
    # Create a dummy input on the same device as the model
    device = next(model.parameters()).device
    dummy_input = torch.randn(1, input_size, device=device)
    
    scaled_model = ScaledNeuralNet(model, scaler, device)
    
    # Make sure model is in evaluation mode
    scaled_model.eval()
    
    # Export the model
    torch.onnx.export(
        scaled_model,                          # model being run
        dummy_input,                           # model input (or a tuple for multiple inputs)
        save_path,                             # where to save the model
        export_params=True,                    # store the trained parameter weights inside the model file
        opset_version=12,                      # the ONNX version to export the model to
        do_constant_folding=True,              # whether to execute constant folding for optimization
        input_names=['input'],                 # the model's input names
        output_names=['output'],               # the model's output names
        dynamic_axes={
            'input': {0: 'batch_size'},        # variable length axes
            'output': {0: 'batch_size'}
        }
    )
    print(f"Model successfully exported to ONNX at {save_path}")

#################################################
################ MAIN FUNCTION ##################
#################################################

def main(data_path, model_save_path=MODEL_NAME):
    
    try:
        with open(data_path, 'r') as f:
            print("Data preview:")
            for i, line in enumerate(f):
                print(line.strip())
                if i >= 5:  # Show first 5 lines
                    break
    except Exception as e:
        print(f"Error reading file: {e}")
    
    device = torch.device("cuda:0" if (not USE_ONLY_CPU and torch.cuda.is_available()) else "cpu")
    print(f"Using device: {device}")

    # Load dataset
    dataset = DubinsDatasetRectangle(data_path, use_trigonometric_features=TRIG_FUNCS)
    print(f"Dataset loaded: {len(dataset)} samples")
    
    # Split dataset
    train_size = int(0.7 * len(dataset))
    val_size = int(0.15 * len(dataset))
    test_size = len(dataset) - train_size - val_size
    
    train_dataset, val_dataset, test_dataset = random_split(
        dataset, [train_size, val_size, test_size]
    )
    
    # Create data loaders
    batch_size = BATCH_SIZE
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size)
    test_loader = DataLoader(test_dataset, batch_size=batch_size)
    
    # Initialize model
    input_size = dataset.get_num_features()  # Number of features
    hidden_size = HIDDEN_SIZE
    output_size = 2  # sin and cos components
    
    model = NeuralNet(input_size, hidden_size, output_size)
    
    # Define loss function and optimizer
    criterion = nn.MSELoss()  # Mean Squared Error loss for regression
    optimizer = optim.Adam(model.parameters(), lr=LEARN_RATE, weight_decay=WEIGHT_DEC)
    
    if DO_TRAINING:
        # Train model       
        print("Starting training...")
        trained_model, history = model_train(
            model, train_loader, val_loader, criterion, optimizer, device, 
            num_epochs=EPOCHS, patience=PATIENCE
        )
        
        # Save the model
        torch.save(trained_model.state_dict(), model_save_path)
        print(f"Model saved to {model_save_path}")
        
        # Plot training results
        plot_training_results(history)
        
        # Evaluate on test set
        print("\nEvaluating on test set...")
        test_loss, _, _ = evaluate_model(trained_model, test_loader, criterion, device)
        
        if EXPORT_TO_ONNX:
            export_to_onnx(trained_model, ONNX_NAME, input_size, dataset.get_scaler())
    
    scaler = dataset.get_scaler()
    inference = ModelInference(model_save_path, scaler, input_size, hidden_size, device)
    
    if EXPORT_TO_ONNX and not DO_TRAINING:
        export_to_onnx(model, ONNX_NAME, input_size, dataset.get_scaler())
    
    testset = DubinsDatasetRectangle(data_path, use_trigonometric_features=TRIG_FUNCS)
    test_loader = DataLoader(testset, batch_size=1024)

    #     # Split dataset
    # train_size = int(0.7 * len(dataset))
    # val_size = int(0.15 * len(dataset))
    # test_size = len(dataset) - train_size - val_size
    
    # train_dataset, val_dataset, test_dataset = random_split(
    #     dataset, [train_size, val_size, test_size]
    # )
    
    # # Create data loaders
    # batch_size = BATCH_SIZE
    # train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)

    test_loss, all_pred, all_tar = evaluate_model(inference.model, test_loader, criterion, device)
    print("Model evaluation complete.")

    # # Example features for prediction
    # # k_max = 1.4
    # # th_i = 1.309 
    # # th_f = -1.11022e-16
    # # alpha_m = -1.5708
    # # alpha_f = 2.87979
    # # th_m = 4.46804
    
    # k_max = 4 
    # th_i = -2.0944
    # th_f = -3.14159 
    # alpha_m = 0.785398 
    # alpha_f = -1.0472 
    # th_m = 0.0174533 
    # man = 11 
    # # len = 3.26627

    # print("th_i=-np.pi/2+np.arctan(0.25)", th_i)
    # print("th_f=-3.0*np.pi/2+np.arctan(0.25)", th_f)
    # print("alpha_m=-np.pi/2", alpha_m)
    # print("alpha_f=2*np.arctan(0.25)-np.pi", alpha_f)
    # print("k_max=np.sqrt(17)/2", k_max)

    # cos_th_i = np.cos(th_i)
    # sin_th_i = np.sin(th_i)
    # cos_th_f = np.cos(th_f)
    # sin_th_f = np.sin(th_f)
    # cos_alpha_m = np.cos(alpha_m)
    # sin_alpha_m = np.sin(alpha_m)
    # cos_alpha_f = np.cos(alpha_f)
    # sin_alpha_f = np.sin(alpha_f)
    
    # example_features = np.array([k_max, sin_th_i, cos_th_i, sin_th_f, cos_th_f, sin_alpha_m, cos_alpha_m, sin_alpha_f, cos_alpha_f])    
    
    # sin_val, cos_val, angle_rad, angle_deg = inference.predict(example_features)
    # error = angle_rad - th_m
    # while error < -np.pi:
    #     error += 2*np.pi
    # while error > np.pi:
    #     error -= 2*np.pi
    # print(f"Error: {error} {np.rad2deg(error)}")

    # print(f"Example features: {example_features}")
    # print(f"Predicted sin(th_m): {sin_val:.4f}, cos(th_m): {cos_val:.4f}")
    # print(f"Predicted angle: {angle_rad:.4f} radians ({angle_deg:.2f} degrees)")
    
    
if __name__ == "__main__":
    main(DATASET_NAME)