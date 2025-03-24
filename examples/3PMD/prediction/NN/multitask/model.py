import torch.nn as nn

class MultiTaskNeuralNet(nn.Module):
    def __init__(self, in_size=5, hid_size=64, class_out_size=19, reg_out_size=1):
        super(MultiTaskNeuralNet, self).__init__()

        # Shared feature extractor
        self.feature_extractor = nn.Sequential(
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
            nn.ReLU()
        )

        # Classification head
        self.classification_head = nn.Sequential(
            nn.Linear(hid_size, class_out_size),
            nn.Softmax(dim=1)
        )

        # Regression head
        self.regression_head = nn.Linear(hid_size, reg_out_size)

        self.initialize_weights()

    def forward(self, x):
        features = self.feature_extractor(x)
        class_out = self.classification_head(features)  
        reg_out = self.regression_head(features)  
        return class_out, reg_out

    def initialize_weights(self):
        for layer in self.feature_extractor:
            if isinstance(layer, nn.Linear):
                nn.init.xavier_normal_(layer.weight)
                nn.init.zeros_(layer.bias)

        for layer in self.classification_head:
            if isinstance(layer, nn.Linear):
                nn.init.xavier_normal_(layer.weight)
                nn.init.zeros_(layer.bias)

        nn.init.xavier_normal_(self.regression_head.weight)
        nn.init.zeros_(self.regression_head.bias)
