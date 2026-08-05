import torch.nn as nn


class NeuralNet(nn.Module):
    def __init__(self, in_size, hid_size, out_size, layers=None):
        super().__init__()

        self.in_size = in_size
        self.hid_size = hid_size
        self.out_size = out_size
        self.layers = layers or []

        if self.layers:
            self.linear_stack = self._build_from_config(self.layers)
        else:
            self.linear_stack = nn.Sequential(
                nn.Linear(in_size, hid_size),
                nn.ReLU(),
                nn.Dropout(0.1),
                nn.Linear(hid_size, hid_size),
                nn.ReLU(),
                nn.Linear(hid_size, out_size),
            )

        self.initialize_weights()

    def _build_from_config(self, layers):
        env = {
            "Linear": nn.Linear,
            "ReLU": nn.ReLU,
            "LeakyReLU": nn.LeakyReLU,
            "Dropout": nn.Dropout,
            "BatchNorm1d": nn.BatchNorm1d,
            "in_size": self.in_size,
            "hid_size": self.hid_size,
            "out_size": self.out_size,
        }
        return nn.Sequential(*[eval(layer, {"__builtins__": {}}, env) for layer in layers])

    def forward(self, x):
        return self.linear_stack(x)

    def initialize_weights(self):
        for layer in self.linear_stack:
            if isinstance(layer, nn.Linear):
                nn.init.xavier_normal_(layer.weight)
                nn.init.zeros_(layer.bias)

    def summarize(self):
        total_params = sum(p.numel() for p in self.parameters() if p.requires_grad)
        print(
            f"NeuralNet(input={self.in_size}, hidden={self.hid_size}, "
            f"output={self.out_size}, trainable_params={total_params})"
        )
        for layer in self.linear_stack:
            print(f"  {layer}")
        return total_params
