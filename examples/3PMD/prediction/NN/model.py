import torch.nn as nn
import torch

class NeuralNet(nn.Module):

    def __init__(self, in_size=5, hid_size=64, out_size=19):
        super(NeuralNet, self).__init__()
        
        self.linear_stack = nn.Sequential(
          nn.Linear(in_size, hid_size),
          nn.ReLU(),
          nn.Linear(hid_size, hid_size),
          nn.ReLU(),
          nn.Linear(hid_size, 128),
          nn.ReLU(),
          nn.Linear(128, 128),
          nn.ReLU(),
          nn.Linear(128, hid_size),
          nn.Dropout(),
          nn.ReLU(),
          nn.Linear(hid_size, out_size),
          nn.Sigmoid()
        )
        
        self.initialize_weights()
        
    def forward(self, x):
        out = self.linear_stack(x)
        return out
    
    def initialize_weights(self):
        for layer in self.linear_stack:
            if isinstance(layer, nn.Linear):
                nn.init.xavier_normal_(layer.weight)
                nn.init.zeros_(layer.bias)