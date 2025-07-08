import torch.nn as nn

class NeuralNet(nn.Module):

    def __init__(self, in_size=5, hid_size=64, out_size=2):  # Changed out_size to 2 for sin and cos
        super(NeuralNet, self).__init__()
        
        self.linear_stack = nn.Sequential(
          nn.Linear(in_size, hid_size // 2),
          nn.ReLU(),
          nn.Dropout(0.2),
          nn.Linear(hid_size // 2, hid_size),
          nn.LeakyReLU(),
          nn.Linear(hid_size, hid_size),
          nn.LeakyReLU(),
          nn.Linear(hid_size, hid_size // 2),
          nn.ReLU(),
          nn.Dropout(0.2),
          nn.Linear(hid_size // 2, out_size)
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


# -- 5 layers 
# 128 - Test MSE 0.3318 MAE 0.1601
# 256 - Test MSE 0.3017 MAE 0.1430
# 512 - Test MSE 0.2766 MAE 0.1252

# -- 4 layers
# 512 - Test MSE 0.3252 MAE 0.1615

# -- 6 layers
# 128 - Test MSE 0.3363 MAE 0.1611
# 512 - Test MSE 0.2718 MAE 0.1244