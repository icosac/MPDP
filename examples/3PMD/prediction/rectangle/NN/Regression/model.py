import torch.nn as nn

class NeuralNet(nn.Module):
    def __init__(self, in_size, hid_size, out_size, layers=[]):
        super(NeuralNet, self).__init__()

        self.in_size = in_size
        self.hid_size = hid_size
        self.out_size = out_size

        if len(layers) > 0:
            self.__build_from_config(layers)
        else:
            print("Using default architecture")
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
            self.print_model()

        self.initialize_weights()


    def __build_from_config(self, layers):
        config = {"layers": layers}

        self.__print_yaml_config(config)

        # Safe evaluation context
        env = {
            "Linear": nn.Linear,
            "ReLU": nn.ReLU,
            "LeakyReLU": nn.LeakyReLU,
            "Dropout": nn.Dropout,
            "in_size": self.in_size,
            "hid_size": self.hid_size,
            "out_size": self.out_size
        }

        self.linear_stack = nn.Sequential(
            *[eval(layer_str, {}, env) for layer_str in config["layers"]]
        )


    def forward(self, x):
        return self.linear_stack(x)

    def initialize_weights(self):
        for layer in self.linear_stack:
            if isinstance(layer, nn.Linear):
                nn.init.xavier_normal_(layer.weight)
                nn.init.zeros_(layer.bias)

    
    def print_model(self):
        print("NeuralNet Architecture: input size {}, hidden size {}, output size {}".format(self.in_size, self.hid_size, self.out_size))
        for layer in self.linear_stack:
            if isinstance(layer, nn.Linear):
                print(f"Linear({layer.in_features}, {layer.out_features})")
            elif isinstance(layer, nn.ReLU):
                print("ReLU()")
            elif isinstance(layer, nn.Dropout):
                print(f"Dropout({layer.p})")
            else:
                print(f"Unknown Layer: {layer}")
        print("For Excel: ", end="")
        for layer in self.linear_stack:
            if isinstance(layer, nn.Linear):
                print(f"L({layer.in_features}, {layer.out_features})", end=",")
            elif isinstance(layer, nn.ReLU):
                print("ReLu()", end=",")
            elif isinstance(layer, nn.Dropout):
                print(f"D({layer.p})", end=",")
        print()


    def __print_yaml_config(self, config):
        print("Layers:")
        for layer in config["layers"]:
            print(f"  - {layer}")
        print("For Excel: ", end="")
        for layer in config["layers"]:
            print_layer = layer.strip()
            if 'hidden_size' in print_layer or 'hid_size' in print_layer:
                print_layer = print_layer.replace("hidden_size", "HS").replace("hid_size", "HS")
            if 'in_size' in print_layer:
                print_layer = print_layer.replace("in_size", "I")
            if 'out_size' in print_layer:
                print_layer = print_layer.replace("out_size", "O")
            if 'Linear' in print_layer:
                print_layer = print_layer.replace("Linear", "L")
            if 'Dropout' in print_layer:
                print_layer = print_layer.replace("Dropout", "D")
            print(print_layer, end=",")
        print()  