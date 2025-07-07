import torch
import os
from model import NeuralNet
import argparse
import numpy as np

PROJECT_PATH = os.path.dirname(os.path.abspath(__file__))
SAVE_PATH = os.path.join(PROJECT_PATH, "models")
SAVE_NAME = os.path.join(SAVE_PATH, "model.pt")

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
device = "cpu"
print("Using device:", device)

def predict_top_n(model, X_input, N=3):
    model.eval()
    X_input = X_input.to(device) 

    with torch.no_grad():
        logits = model(X_input)
        # For multi-label, use raw logits and top-k
        top_n_vals, top_n_classes = torch.topk(logits, N, dim=1)

    return top_n_classes.cpu().numpy(), top_n_vals.cpu().numpy()

def preprocess_input(data):
    """
    Preprocess input data by computing cosine and sine for angle features.
    """
    angles = ['theta_i', 'theta_f', 'alpha_m', 'alpha_f']
    processed_data = []
    for i, value in enumerate(data):
        if angles[i % len(angles)] in angles:
            angle_rad = np.radians(value)
            processed_data.append(np.cos(angle_rad))
            processed_data.append(np.sin(angle_rad))
        else:
            processed_data.append(value)
    return processed_data

def preprocess_data(data):
    """
    Preprocess data by computing cosine and sine for angle features.
    """
    angles = ['theta_i', 'theta_f', 'alpha_m', 'alpha_f']
    for angle in angles:
        data[f'{angle}_cos'] = np.cos(np.radians(data[angle]))
        data[f'{angle}_sin'] = np.sin(np.radians(data[angle]))
    return data

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Inference script for the Neural Network model')
    parser.add_argument('--input', type=str, help='Input data for prediction', required=True)
    args = parser.parse_args()
    
    model = NeuralNet()  
    model.load_state_dict(torch.load(SAVE_NAME, map_location=device, weights_only=True))
    model.to(device)

    if args.input.endswith(".csv"):
        import pandas as pd
        data = pd.read_csv(args.input, sep='\s+')
        
        original_features = ['kmax', 'theta_i', 'theta_f', 'alpha_m', 'alpha_f', 'id_man_comb']
        data = data[original_features]

        data = preprocess_data(data)

        # Update feature list to include cosine and sine of angles
        features = ['kmax', 'theta_i_cos', 'theta_i_sin', 'theta_f_cos', 'theta_f_sin',
                    'alpha_m_cos', 'alpha_m_sin', 'alpha_f_cos', 'alpha_f_sin']
        target = 'id_man_comb'

        # Preprocess data
        X = data[features].values
        y = data[target].values

        tot = 0 
        correct = 0

        for entry, label in zip(X, y):
            data_tensor = torch.tensor([entry], dtype=torch.float32)
            top_classes, top_vals = predict_top_n(model, data_tensor, N=3)
            # Parse label as set of valid classes (multi-label)
            if isinstance(label, str):
                valid_classes = set(int(i) for i in label.split(',') if i.strip() != '')
            else:
                valid_classes = {int(label)}
            # Check if any predicted class is in the set of valid classes
            if any(cls in valid_classes for cls in top_classes[0]):
                correct += 1
            tot += 1

        print("Top-3 accuracy (any predicted class in label set):", correct / tot * 100.0)

    else:
        # Transform string to list of floats
        args.input = args.input.split(',')
        args.input = [float(i) for i in args.input]

        # Preprocess input to compute cosine and sine for angles
        args.input = preprocess_input(args.input)
        
        data = torch.tensor([args.input], dtype=torch.float32)
        
    top_classes, top_probs = predict_top_n(model, data, N=5)

    print("Top-N Predicted Classes:", top_classes)
    print("Top-N Probabilities:", top_probs)
    print("Top-N Probabilities:", top_probs)
