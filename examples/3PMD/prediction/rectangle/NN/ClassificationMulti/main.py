#################################################
############### CLASSIFICATION ##################
#################################################

import torch
import torch.onnx
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, random_split
import matplotlib.pyplot as plt
from model import NeuralNet
from dataset import DubinsDatasetRectangle
import numpy as np
from train import model_train, evaluate_model

import os
from pathlib import Path

BATCH_SIZE  = 32
EPOCHS      = 15
PATIENCE    = 10
LEARN_RATE  = 0.0001
WEIGHT_DEC  = 1e-5
HIDDEN_SIZE = 1024

TRIG_FUNCS  = True  # Use trigonometric features

DO_TRAINING    = True
EXPORT_TO_ONNX = False

THIS_FILE_PATH = os.path.abspath(__file__)
PROJECT_PATH   = Path(THIS_FILE_PATH).parent
DATASET_PATH   = os.path.join(PROJECT_PATH.parent.parent, "datasets")
MODELS_PATH    = os.path.join(PROJECT_PATH, "models")
PLOT_PATH      = os.path.join(PROJECT_PATH, "plots")

DATASET_NAME   = os.path.join(DATASET_PATH, "/Users/enrico/Projects/mpdp/small_rec.csv")
MODEL_NAME     = os.path.join(MODELS_PATH, 'model.pt')
ONNX_NAME      = os.path.join(MODELS_PATH, 'model.onnx')

os.makedirs(MODELS_PATH, exist_ok=True)
os.makedirs(PLOT_PATH,   exist_ok=True)

#################################################
################ PLOT FUNCS #####################
#################################################

def plot_training_results(history):
    """
    Plots training and validation metrics
    
    Args:
        history: Dictionary containing training history
    """
    plt.figure(figsize=(12, 5))
    
    # Plot loss
    plt.subplot(1, 2, 1)
    plt.plot(history['train_losses'], label='Train Loss')
    plt.plot(history['val_losses'], label='Validation Loss')
    plt.title('Loss over Epochs')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    
    # Plot accuracy
    plt.subplot(1, 2, 2)
    plt.plot(history['train_accs'], label='Train Accuracy')
    plt.plot(history['val_accs'], label='Validation Accuracy')
    plt.title('Accuracy over Epochs')
    plt.xlabel('Epoch')
    plt.ylabel('Accuracy')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_PATH, 'training_history.png'))

#################################################
################ INFERENCE CLASS ################
#################################################

class ModelInference:
    def __init__(self, model_path, scaler, input_size, hidden_size, num_classes, device='cpu', inverse_mapping=None):
        """
        Class for inference with a trained model.
        
        Args:
            model_path: Path to the saved model
            scaler: StandardScaler fit on training data
            input_size: Number of input features
            hidden_size: Number of hidden units
            num_classes: Number of output classes
            device: Device to use for computation
            inverse_mapping: Dictionary to map model output indices back to original class labels
        """
        self.device = device
        self.scaler = scaler
        self.inverse_mapping = inverse_mapping
        
        # Initialize model
        self.model = NeuralNet(input_size, hidden_size, num_classes)
        self.model.load_state_dict(torch.load(model_path, map_location=device))
        self.model.to(device)
        self.model.eval()
        
    def predict(self, features, top_k=4):
        """
        Predicts the class of a sample.
        
        Args:
            features: Features as numpy array shape (5,) or (n, 5)
            top_k: Number of top classes to return (for multi-label problems)
        
        Returns:
            Tuple of (predicted_class, original_class_label) or arrays of these
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
            # Get top-k results for multi-label
            topk_vals, topk_indices = torch.topk(outputs, top_k, dim=1)
        model_output = topk_indices.cpu().numpy()
        
        # Map back to original labels if mapping exists
        if self.inverse_mapping:
            original_labels = np.array([[self.inverse_mapping[idx] for idx in row] for row in model_output])
        else:
            original_labels = model_output
        
        if single_sample:
            return model_output[0], original_labels[0]
        else:
            return model_output, original_labels
    


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
        dummy_input,                    # model input (or a tuple for multiple inputs)
        save_path,                      # where to save the model
        export_params=True,             # store the trained parameter weights inside the model file
        opset_version=12,               # the ONNX version to export the model to
        do_constant_folding=True,       # whether to execute constant folding for optimization
        input_names=['input'],          # the model's input names
        output_names=['output'],        # the model's output names
        dynamic_axes={
            'input': {0: 'batch_size'},  # variable length axes
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
    
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # Load dataset
    dataset = DubinsDatasetRectangle(data_path, use_trigonometric_features=TRIG_FUNCS)
    print(f"Dataset loaded: {len(dataset)} samples, {dataset.num_classes} classes")
    
    # Split dataset
    train_size = int(0.7 * len(dataset))
    val_size = int(0.15 * len(dataset))
    test_size = len(dataset) - train_size - val_size
    
    train_dataset, val_dataset, test_dataset = random_split(
        dataset, [train_size, val_size, test_size]
    )
    
    # Create data loaders
    train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=BATCH_SIZE)
    test_loader = DataLoader(test_dataset, batch_size=BATCH_SIZE)
    
    # Initialize model
    input_size = dataset.get_num_features()
    num_classes = dataset.num_classes

    print(f"Input size: {input_size}, Hidden size: {HIDDEN_SIZE}, Num classes: {num_classes}")
    
    model = NeuralNet(input_size, HIDDEN_SIZE, num_classes)
    
    # Define loss function and optimizer
    criterion = nn.BCEWithLogitsLoss()
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
        test_acc, _, _ = evaluate_model(trained_model, test_loader, criterion, device)
    
    scaler = dataset.get_scaler()
    label_mapping, inverse_mapping = dataset.get_label_mapping()
    inference = ModelInference(model_save_path, scaler, input_size, HIDDEN_SIZE, num_classes, device, inverse_mapping)
    
    # Class of this sample is 12
    if TRIG_FUNCS:
        from math import sin, cos
        example_features = np.array([1, 0.1, 0, 0.1, 0.5, sin(-0.2618), cos(-0.2618), sin(1.5708), cos(1.5708)])  
    else:
        example_features = np.array([1, 0.1, 0, 0.1, 0.5, -0.2618, 1.5708])

    # Predict
    predicted_classes, original_classes = inference.predict(example_features)
    print(f"Example features: {example_features}")
    print(f"Top-4 predicted class indices: {predicted_classes}")
    print(f"Top-4 original class IDs (in your data): {original_classes}")
    print("Note: The label for each sample is a set of valid classes. Prediction is considered correct if any predicted class is in the set.")
    
    if EXPORT_TO_ONNX:
        export_to_onnx(model, ONNX_NAME, input_size, dataset.get_scaler())


if __name__ == "__main__":
    main(DATASET_NAME)
    main(DATASET_NAME)
