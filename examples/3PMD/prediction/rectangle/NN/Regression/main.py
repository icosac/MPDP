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
import time
import yaml
import argparse

RANDOM_SEED = 42

BATCH_SIZE  = 1
EPOCHS      = 1
PATIENCE    = 1
WEIGHT_DEC  = 1
LEARN_RATE  = 1

LAYERS      = []  
HIDDEN_SIZE = 1

TRIG_FUNCS  = True  # Use trigonometric features

DO_TRAINING       = True
USE_ONLY_CPU      = False
EXPORT_TO_ONNX    = True

THIS_FILE_PATH    = os.path.abspath(__file__)
PROJECT_PATH      = Path(THIS_FILE_PATH).parent
DATASET_PATH      = os.path.join(PROJECT_PATH.parent.parent, "datasets")
OUTPUT_MODEL_PATH = os.path.join(PROJECT_PATH, "models")
PLOT_PATH         = os.path.join(PROJECT_PATH, "plots")

SLURM_JOB_ID      = os.environ.get('SLURM_JOB_ID', '')
SLURM_ID_STR      = f"_{SLURM_JOB_ID}" if SLURM_JOB_ID else ""

DATASET_NAME      = os.path.join("/Users/enrico/Projects/mpdp/small_rect.csv")
MODEL_NAME        = os.path.join(OUTPUT_MODEL_PATH, 'model_regression{}.pt'.format(SLURM_ID_STR))
ONNX_NAME         = os.path.join(OUTPUT_MODEL_PATH, 'model_regression{}.onnx'.format(SLURM_ID_STR))


#################################################
################ PLOT FUNCS #####################
#################################################

def plot_training_results(history):
    """
    Plots training and validation metrics
    
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
    plt.savefig(os.path.join(PLOT_PATH, 'training_history{}.png'.format(SLURM_ID_STR)))

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
        self.model = NeuralNet(input_size, hidden_size, out_size=2, layers=LAYERS)
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
    train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=BATCH_SIZE)
    test_loader = DataLoader(test_dataset, batch_size=BATCH_SIZE)
    
    input_size = dataset.get_num_features()  # Number of features
    hidden_size = HIDDEN_SIZE
    output_size = 2  # sin and cos components
    
    if DO_TRAINING:
        # Initialize model
        model = NeuralNet(input_size, hidden_size, output_size, layers=LAYERS)
        
        # Define loss function and optimizer
        criterion = nn.MSELoss()  # Mean Squared Error loss for regression
        optimizer = optim.Adam(model.parameters(), lr=LEARN_RATE, weight_decay=WEIGHT_DEC)

        training_time = time.time()
        # Train model       
        print("Starting training...")
        trained_model, history = model_train(
            model, train_loader, val_loader, criterion, optimizer, device, 
            num_epochs=EPOCHS, patience=PATIENCE
        )

        print(f"Training completed in {time.time() - training_time:.4f} seconds")
        training_time = time.time()

        # Save the model
        torch.save(trained_model.state_dict(), model_save_path)
        if EXPORT_TO_ONNX:
            export_to_onnx(trained_model, ONNX_NAME, input_size, dataset.get_scaler())
        print(f"Model saved to {model_save_path} in {time.time() - training_time:.4f} seconds")
        
        # Plot training results
        plot_training_results(history)
        
        # Evaluate on test set
        print("\nEvaluating on test set...")
        training_time = time.time()
        test_loss, _, _ = evaluate_model(trained_model, test_loader, criterion, device, plot_path=PLOT_PATH, slurm_id_str=SLURM_ID_STR, eval_test=True)
        print(f"Model evaluation complete in {time.time() - training_time:.4f} seconds")
    
    else:
        print("Running inference...")
        scaler = dataset.get_scaler()
        inference_model = ModelInference(
            model_path=model_save_path, 
            scaler=scaler, 
            input_size=input_size, 
            hidden_size=hidden_size, 
            device=device
        )

        if EXPORT_TO_ONNX:
            print("Exporting to ONNX")
            export_to_onnx(model, ONNX_NAME, input_size, dataset.get_scaler())
        
        sin, cos, rad, deg = inference_model.predict(np.array([2, 1, 0.25, 0.75, np.sin(np.pi/2.0), np.cos(np.pi/2.0), np.sin(-np.pi/2.0), np.cos(-np.pi/2.0)]))

        print(sin, cos)
        print(2*np.pi+rad, 360+deg)

        # testset = DubinsDatasetRectangle(data_path, use_trigonometric_features=TRIG_FUNCS)
        # test_loader = DataLoader(testset, batch_size=1024)

        # test_loss, _, _ = evaluate_model(inference_model.model, test_loader, criterion, device, plot_path=PLOT_PATH, slurm_id_str=SLURM_ID_STR, eval_test=True)
        # print("Model evaluation complete.")

    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run regression model for 3 Point Dubins path prediction.")

    parser.add_argument('--yaml-config', type=str, default=None, help='Path to YAML configuration file. If additional arguments are provided, they will override the YAML config.')

    # Set all arguments to default=None except yaml-config
    parser.add_argument('--random-seed', type=int, default=None, help='Random seed for reproducibility')
    parser.add_argument('--batch-size', type=int, default=None, help='Batch size for training and evaluation')
    parser.add_argument('--epochs', type=int, default=None, help='Number of epochs for training')
    parser.add_argument('--patience', type=int, default=None, help='Early stopping patience')
    parser.add_argument('--weight-decay', type=float, default=None, help='Weight decay for the optimizer')
    parser.add_argument('--learn-rate', type=float, default=None, help='Learning rate for the optimizer')
    parser.add_argument('--hidden-size', type=int, default=None, help='Number of hidden units in the model')
    parser.add_argument('--trig-funcs', type=lambda x: (str(x).lower() == 'true'), default=None, help='Enable trigonometric features (sin, cos)')
    parser.add_argument('--do-training', type=lambda x: (str(x).lower() == 'true'), default=None, help='Enable training mode')
    parser.add_argument('--use-only-cpu', type=lambda x: (str(x).lower() == 'true'), default=None, help='Use only CPU for training and evaluation')
    parser.add_argument('--export-to-onnx', type=lambda x: (str(x).lower() == 'true'), default=None, help='Export model to ONNX format')
    parser.add_argument('--dataset', type=str, default=None, help='Path to the dataset CSV file')
    parser.add_argument('--output-model-path', type=str, default=None, help='Path to save the trained model')
    parser.add_argument('--output-model-name', type=str, default=None, help='Name of the output model file')
    parser.add_argument('--onnx-name', type=str, default=None, help='Name of the output ONNX file')
    parser.add_argument('--plot-path', type=str, default=None, help='Path to save training plots')

    args = parser.parse_args()

    # Load YAML config if provided
    config = {}
    if args.yaml_config:
        with open(args.yaml_config, 'r') as f:
            config = yaml.safe_load(f)
            print(f"Loaded configuration from {args.yaml_config}")
            print(config)

    # Helper function to resolve value priority
    def resolve_arg(arg_val, config_key, default_val):
        if arg_val is not None:
            return arg_val
        elif config_key in config:
            if config_key in ['BATCH_SIZE', 'EPOCHS', 'PATIENCE', 'HIDDEN_SIZE']:
                return int(config[config_key])
            elif config_key in ['LEARN_RATE', 'WEIGHT_DEC']:
                return float(config[config_key])
            elif config_key in ['TRIG_FUNCS', 'DO_TRAINING', 'USE_ONLY_CPU', 'EXPORT_TO_ONNX']:
                return bool(config[config_key])
            else:
                return config[config_key]
        else:
            return default_val
           
    # Set all parameters with correct priority
    RANDOM_SEED = resolve_arg(args.random_seed, 'RANDOM_SEED', RANDOM_SEED)
    BATCH_SIZE = resolve_arg(args.batch_size, 'BATCH_SIZE', BATCH_SIZE)
    EPOCHS = resolve_arg(args.epochs, 'EPOCHS', EPOCHS)
    PATIENCE = resolve_arg(args.patience, 'PATIENCE', PATIENCE)
    WEIGHT_DEC = resolve_arg(args.weight_decay, 'WEIGHT_DEC', WEIGHT_DEC)
    LEARN_RATE = resolve_arg(args.learn_rate, 'LEARN_RATE', LEARN_RATE)
    LAYERS = resolve_arg(None, 'LAYERS', LAYERS)
    HIDDEN_SIZE = resolve_arg(args.hidden_size, 'HIDDEN_SIZE', HIDDEN_SIZE)
    TRIG_FUNCS = resolve_arg(args.trig_funcs, 'TRIG_FUNCS', TRIG_FUNCS)
    DO_TRAINING = resolve_arg(args.do_training, 'DO_TRAINING', DO_TRAINING)
    USE_ONLY_CPU = resolve_arg(args.use_only_cpu, 'USE_ONLY_CPU', USE_ONLY_CPU)
    EXPORT_TO_ONNX = resolve_arg(args.export_to_onnx, 'EXPORT_TO_ONNX', EXPORT_TO_ONNX)
    DATASET_NAME = resolve_arg(args.dataset, 'DATASET', DATASET_NAME)
    OUTPUT_MODEL_PATH = resolve_arg(args.output_model_path, 'OUTPUT_MODEL_PATH', OUTPUT_MODEL_PATH)
    PLOT_PATH = resolve_arg(args.plot_path, 'PLOT_PATH', PLOT_PATH)
    MODEL_NAME = resolve_arg(args.output_model_name, 'OUTPUT_MODEL_NAME', os.path.join(OUTPUT_MODEL_PATH, 'model_regression{}.pt'.format(SLURM_ID_STR)))
    ONNX_NAME = resolve_arg(args.onnx_name, 'ONNX_NAME', os.path.join(OUTPUT_MODEL_PATH, 'model_regression{}.onnx'.format(SLURM_ID_STR)))

    if not os.path.exists(DATASET_NAME):
        raise FileNotFoundError(f"Dataset file not found: {DATASET_NAME}")

    print(f"Running main with dataset: {DATASET_NAME}")
    print(f"ONNX model will be saved to: {ONNX_NAME}")
    print(f"Plots will be saved to: {PLOT_PATH}")
    print(f"SLURM ID: {SLURM_JOB_ID if SLURM_JOB_ID else 'Not running on SLURM'}")
    print(f"Using device: {torch.device('cuda:0' if (not USE_ONLY_CPU and torch.cuda.is_available()) else 'cpu')}")
    print(f"Trigonometric features: {'Enabled' if TRIG_FUNCS else 'Disabled'}")
    print(f"Training enabled: {DO_TRAINING}")
    print(f"Exporting to ONNX: {'Enabled' if EXPORT_TO_ONNX else 'Disabled'}")
    print(f"Batch size: {BATCH_SIZE}, Epochs: {EPOCHS}, Patience: {PATIENCE}")
    print(f"Learning rate: {LEARN_RATE}, Weight decay: {WEIGHT_DEC}, Hidden size: {HIDDEN_SIZE}")
    print(f"Dataset path: {DATASET_NAME}")
    print(f"Model save path: {MODEL_NAME}")
    print(f"ONNX save path: {ONNX_NAME}")
    print(f"Plot path: {PLOT_PATH}")
    print(f"Random seed: {RANDOM_SEED}")

    # Set random seed for reproducibility
    torch.manual_seed(RANDOM_SEED)
    np.random.seed(RANDOM_SEED)

    os.makedirs(OUTPUT_MODEL_PATH, exist_ok=True)
    os.makedirs(PLOT_PATH, exist_ok=True)

    main(DATASET_NAME, model_save_path=MODEL_NAME)   