import torch.nn as nn

class NeuralNet(nn.Module):

    def __init__(self, in_size=9, hid_size=64, out_size=19):  # Revert to original: out_size is total number of classes
        super(NeuralNet, self).__init__()
        self.linear_stack = nn.Sequential(
            nn.Linear(in_size, hid_size // 2),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(hid_size // 2, hid_size),
            nn.ReLU(),
            nn.Linear(hid_size, hid_size),
            nn.ReLU(),
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
