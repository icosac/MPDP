import torch
import os
from model import NeuralNet
import argparse
from sklearn.preprocessing import OneHotEncoder

import numpy as np
np.set_printoptions(threshold=np.inf)

PROJECT_PATH = os.path.dirname(os.path.abspath(__file__))
SAVE_PATH = os.path.join(PROJECT_PATH, "models")
SAVE_NAME = os.path.join(SAVE_PATH, "model.pt")

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Using device:", device)

def predict_entry(model, X_test, N=3):
    model.eval()
    with torch.no_grad():
        X_test = X_test.to(device)
        y_prob, y_pred = torch.topk(torch.softmax(model(X_test), dim=0), N, dim=0)
    return y_pred.cpu().numpy(), y_prob.cpu().numpy()

def predict_batch(model, X_test, N):
    model.eval()
    with torch.no_grad():
        X_test = X_test.to(device)
        y_prob, y_pred = torch.topk(torch.softmax(model(X_test), dim=1), N, dim=1)
        print(y_prob)
    return y_pred.cpu().numpy(), y_prob.cpu().numpy()

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
        
        features = ['kmax', 'theta_i', 'theta_f', 'alpha_m', 'alpha_f']
        target = 'id_man_comb'

        X = data[features].values
        y = data[target].values

        encoder = OneHotEncoder(sparse_output=False)  # Use dense array output
        y = encoder.fit_transform(y.reshape(-1, 1))  # Convert y to numpy and reshape
        
        # Preprocess data
        X = torch.tensor(X, dtype=torch.float32)
        y = torch.tensor(y, dtype=torch.float32)

        top_classes, top_probs = predict_batch(model, X, N=5)

        y = torch.argmax(y, dim=1).numpy()
        
        tot = 0
        count = 0

        for top_class, label in zip(top_classes, y):
            if label in top_class:
                count += 1
            tot += 1

        print("Accuracy:", count / tot * 100.0)

    else:
        # Transform string to list of floats
        args.input = args.input.split(',')
        args.input = [float(i) for i in args.input]

        # Preprocess input to compute cosine and sine for angles
        
        data = torch.tensor([args.input], dtype=torch.float32)
        
        top_classes, top_probs = predict_entry(model, data, N=5)

        print("Top-N Predicted Classes:", top_classes)
        print("Top-N Probabilities:", top_probs)
