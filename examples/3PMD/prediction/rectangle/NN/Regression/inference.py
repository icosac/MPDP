import torch
import os
from model import NeuralNet
import argparse

PROJECT_PATH = os.path.dirname(os.path.abspath(__file__))
SAVE_PATH = os.path.join(PROJECT_PATH, "models")
SAVE_NAME = os.path.join(SAVE_PATH, "model.pt")

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Using device:", device)
device="cpu"

def predict_top_n(model, X_input, N=3):
    model.eval()
    X_input = X_input.to(device) 

    with torch.no_grad():
        logits = model(X_input)  
        probs = torch.softmax(logits, dim=1) 

    top_n_probs, top_n_classes = torch.topk(probs, N, dim=1)  

    return top_n_classes.cpu().numpy(), top_n_probs.cpu().numpy()

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Inference script for the Neural Network model')
    parser.add_argument('--input', type=str, help='Input data for prediction', required=True)
    args = parser.parse_args()
    
    model = NeuralNet()  
    model.load_state_dict(torch.load(SAVE_NAME, map_location=device, weights_only=True))
    model.to(device)

    # transform string to list of floats
    args.input = args.input.split(',')
    args.input = [float(i) for i in args.input]
    
    input = torch.tensor([args.input], dtype=torch.float32)
        
    top_classes, top_probs = predict_top_n(model, input, N=3)

    print("Top-N Predicted Classes:", top_classes)
    print("Top-N Probabilities:", top_probs)