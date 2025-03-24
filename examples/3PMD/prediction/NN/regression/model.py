import torch.nn as nn
import torch

class NeuralNet(nn.Module):
    def __init__(self, in_size=8, hid_size=2048, out_size=2):  # Output is now (sin_th_m, cos_th_m)
        super(NeuralNet, self).__init__()
        
        self.linear_stack = nn.Sequential(
            nn.Linear(in_size, hid_size),
            nn.ReLU(),
            nn.Linear(hid_size, 128),
            nn.ReLU(),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, out_size)  # Outputs two values for angle representation
        )
        
        self.initialize_weights()

    def forward(self, x):
        return self.linear_stack(x)

    def initialize_weights(self):
        for layer in self.linear_stack:
            if isinstance(layer, nn.Linear):
                nn.init.xavier_normal_(layer.weight)
                nn.init.zeros_(layer.bias)
