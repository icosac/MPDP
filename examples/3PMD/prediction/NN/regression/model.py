# import torch.nn as nn

# class NeuralNet(nn.Module):

#     def __init__(self, in_size=5, hid_size=24, out_size=1):
#         super(NeuralNet, self).__init__()
        
#         self.linear_stack = nn.Sequential(
#           nn.Linear(in_size, hid_size),
#           nn.ReLU(),
#           nn.Linear(hid_size, 12),
#           nn.ReLU(),
#           nn.Linear(12, 6),
#           nn.ReLU(),
#           nn.Linear(6, out_size)
#         )
        
#         self.initialize_weights()
        
#     def forward(self, x):
#         out = self.linear_stack(x)
#         return out
    
#     def initialize_weights(self):
#         for layer in self.linear_stack:
#             if isinstance(layer, nn.Linear):
#                 nn.init.xavier_normal_(layer.weight)
#                 nn.init.zeros_(layer.bias)




import torch.nn as nn

class NeuralNet(nn.Module):

    def __init__(self, in_size=5, hid_size=24, out_size=1):
        super(NeuralNet, self).__init__()
        
        self.linear_stack = nn.Sequential(
          nn.Linear(in_size, 100),
          nn.ReLU(),
          nn.Linear(100, 10000),
          nn.ReLU(),
          nn.Dropout(0.2),
          nn.Linear(10000, 100000),
          nn.ReLU(),
          nn.Dropout(0.6),
          nn.Linear(100000, 1000),
          nn.ReLU(),
          nn.Linear(1000, out_size)
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