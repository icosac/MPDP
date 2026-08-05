import argparse
import os
import time
from pathlib import Path

import joblib
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import yaml
from torch.utils.data import DataLoader

from dataset import RSTriangularDataset
from model import NeuralNet
from train import evaluate_model, model_train, plot_training_results


THIS_FILE_PATH = Path(__file__).resolve()
PROJECT_PATH = THIS_FILE_PATH.parent
DATASET_PATH = PROJECT_PATH.parent.parent / "datasets"

DEFAULT_TRAIN_DATASET = DATASET_PATH / "DS_RS_triangular_training.csv"
DEFAULT_VALIDATION_DATASET = DATASET_PATH / "DS_RS_triangular_validation.csv"
DEFAULT_OUTPUT_MODEL_PATH = PROJECT_PATH / "models"
DEFAULT_PLOT_PATH = PROJECT_PATH / "plots"

DEFAULT_CONFIG = {
    "RANDOM_SEED": 42,
    "BATCH_SIZE": 128,
    "EPOCHS": 80,
    "PATIENCE": 10,
    "WEIGHT_DEC": 1e-5,
    "LEARN_RATE": 1e-3,
    "HIDDEN_SIZE": 128,
    "TRIG_FUNCS": True,
    "USE_ONLY_CPU": False,
    "EXPORT_TO_ONNX": True,
    "TRAIN_DATASET": str(DEFAULT_TRAIN_DATASET),
    "VALIDATION_DATASET": str(DEFAULT_VALIDATION_DATASET),
    "OUTPUT_MODEL_PATH": str(DEFAULT_OUTPUT_MODEL_PATH),
    "OUTPUT_MODEL_NAME": str(DEFAULT_OUTPUT_MODEL_PATH / "model_classification.pt"),
    "ONNX_NAME": str(DEFAULT_OUTPUT_MODEL_PATH / "model_classification.onnx"),
    "PLOT_PATH": str(DEFAULT_PLOT_PATH),
    "LAYERS": [
        "Linear(in_size, hid_size)",
        "ReLU()",
        "Dropout(0.1)",
        "Linear(hid_size, hid_size)",
        "ReLU()",
        "Linear(hid_size, out_size)",
    ],
}


class ScaledNeuralNet(nn.Module):
    def __init__(self, model, scaler, use_trigonometric_features=True):
        super().__init__()
        self.model = model
        self.use_trigonometric_features = use_trigonometric_features
        self.register_buffer("mean", torch.tensor(scaler.mean_, dtype=torch.float32))
        self.register_buffer("scale", torch.tensor(scaler.scale_, dtype=torch.float32))

    def forward(self, raw_x):
        if self.use_trigonometric_features:
            thi = raw_x[:, 0:1]
            thf = raw_x[:, 1:2]
            kmax = raw_x[:, 2:3]
            x = torch.cat((torch.sin(thi), torch.cos(thi), torch.sin(thf), torch.cos(thf), kmax), dim=1)
        else:
            x = raw_x

        x = (x - self.mean) / self.scale
        return self.model(x)


def export_to_onnx(model, save_path, scaler, use_trigonometric_features):
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    device = next(model.parameters()).device
    scaled_model = ScaledNeuralNet(model, scaler, use_trigonometric_features).to(device)
    scaled_model.eval()
    dummy_input = torch.randn(1, 3, device=device)

    torch.onnx.export(
        scaled_model,
        dummy_input,
        save_path,
        export_params=True,
        opset_version=18,
        do_constant_folding=True,
        input_names=["input"],
        output_names=["output"],
        dynamic_axes={"input": {0: "batch_size"}, "output": {0: "batch_size"}},
    )
    print(f"Model exported to ONNX at {save_path}")


def train_from_config(config):
    torch.manual_seed(config["RANDOM_SEED"])
    np.random.seed(config["RANDOM_SEED"])

    device = torch.device("cuda:0" if (not config["USE_ONLY_CPU"] and torch.cuda.is_available()) else "cpu")
    print(f"Using device: {device}")

    train_dataset = RSTriangularDataset(
        config["TRAIN_DATASET"],
        use_trigonometric_features=config["TRIG_FUNCS"],
        fit_scaler=True,
    )
    val_dataset = RSTriangularDataset(
        config["VALIDATION_DATASET"],
        scaler=train_dataset.get_scaler(),
        label_mapping=train_dataset.get_label_mapping()[0],
        use_trigonometric_features=config["TRIG_FUNCS"],
    )

    train_loader = DataLoader(train_dataset, batch_size=config["BATCH_SIZE"], shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=config["BATCH_SIZE"], shuffle=False)

    input_size = train_dataset.get_num_features()
    num_classes = train_dataset.num_classes
    print(f"Train samples: {len(train_dataset)}, validation samples: {len(val_dataset)}")
    print(f"Features: {train_dataset.get_feature_names()}")
    print(f"Input size: {input_size}, hidden size: {config['HIDDEN_SIZE']}, classes: {num_classes}")
    print(f"Label mapping: {train_dataset.get_label_mapping()[0]}")

    model = NeuralNet(input_size, config["HIDDEN_SIZE"], num_classes, layers=config["LAYERS"])
    model.summarize()

    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=config["LEARN_RATE"], weight_decay=config["WEIGHT_DEC"])

    start = time.time()
    trained_model, history = model_train(
        model,
        train_loader,
        val_loader,
        criterion,
        optimizer,
        device,
        num_epochs=config["EPOCHS"],
        patience=config["PATIENCE"],
    )
    print(f"Training completed in {time.time() - start:.4f} seconds")

    os.makedirs(config["OUTPUT_MODEL_PATH"], exist_ok=True)
    os.makedirs(config["PLOT_PATH"], exist_ok=True)
    os.makedirs(Path(config["OUTPUT_MODEL_NAME"]).parent, exist_ok=True)

    checkpoint = {
        "model_state_dict": trained_model.state_dict(),
        "input_size": input_size,
        "hidden_size": config["HIDDEN_SIZE"],
        "num_classes": num_classes,
        "layers": config["LAYERS"],
        "use_trigonometric_features": config["TRIG_FUNCS"],
        "feature_names": train_dataset.get_feature_names(),
        "raw_feature_columns": ["thi", "thf", "kmax"],
        "label_mapping": train_dataset.get_label_mapping()[0],
        "inverse_mapping": train_dataset.get_label_mapping()[1],
    }
    torch.save(checkpoint, config["OUTPUT_MODEL_NAME"])
    joblib.dump(train_dataset.get_scaler(), str(Path(config["OUTPUT_MODEL_NAME"]).with_suffix(".scaler.joblib")))
    print(f"Model saved to {config['OUTPUT_MODEL_NAME']}")

    if config["EXPORT_TO_ONNX"]:
        export_to_onnx(
            trained_model,
            config["ONNX_NAME"],
            train_dataset.get_scaler(),
            config["TRIG_FUNCS"],
        )

    plot_training_results(history, config["PLOT_PATH"])
    evaluate_model(
        trained_model,
        val_loader,
        criterion,
        device,
        num_classes,
        plot_path=config["PLOT_PATH"],
        split_name="validation",
    )


def parse_args():
    parser = argparse.ArgumentParser(description="Train a point-to-point RS maneuver classifier.")
    parser.add_argument("--yaml-config", type=str, default=None, help="Path to YAML configuration file.")
    parser.add_argument("--train-dataset", type=str, default=None, help="Path to training CSV.")
    parser.add_argument("--validation-dataset", type=str, default=None, help="Path to validation CSV.")
    parser.add_argument("--batch-size", type=int, default=None)
    parser.add_argument("--epochs", type=int, default=None)
    parser.add_argument("--patience", type=int, default=None)
    parser.add_argument("--hidden-size", type=int, default=None)
    parser.add_argument("--learn-rate", type=float, default=None)
    parser.add_argument("--weight-decay", type=float, default=None)
    parser.add_argument("--trig-funcs", type=str_to_bool, default=None)
    parser.add_argument("--use-only-cpu", type=str_to_bool, default=None)
    parser.add_argument("--export-to-onnx", type=str_to_bool, default=None)
    parser.add_argument("--output-model-path", type=str, default=None)
    parser.add_argument("--output-model-name", type=str, default=None)
    parser.add_argument("--onnx-name", type=str, default=None)
    parser.add_argument("--plot-path", type=str, default=None)
    return parser.parse_args()


def str_to_bool(value):
    if isinstance(value, bool):
        return value
    value = value.lower()
    if value in {"true", "1", "yes", "y"}:
        return True
    if value in {"false", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Expected a boolean value, got {value}")


def load_config(args):
    config = DEFAULT_CONFIG.copy()
    if args.yaml_config:
        with open(args.yaml_config, "r", encoding="utf-8") as stream:
            yaml_config = yaml.safe_load(stream) or {}
        config.update(yaml_config)

    overrides = {
        "TRAIN_DATASET": args.train_dataset,
        "VALIDATION_DATASET": args.validation_dataset,
        "BATCH_SIZE": args.batch_size,
        "EPOCHS": args.epochs,
        "PATIENCE": args.patience,
        "HIDDEN_SIZE": args.hidden_size,
        "LEARN_RATE": args.learn_rate,
        "WEIGHT_DEC": args.weight_decay,
        "TRIG_FUNCS": args.trig_funcs,
        "USE_ONLY_CPU": args.use_only_cpu,
        "EXPORT_TO_ONNX": args.export_to_onnx,
        "OUTPUT_MODEL_PATH": args.output_model_path,
        "OUTPUT_MODEL_NAME": args.output_model_name,
        "ONNX_NAME": args.onnx_name,
        "PLOT_PATH": args.plot_path,
    }
    for key, value in overrides.items():
        if value is not None:
            config[key] = value

    config["OUTPUT_MODEL_PATH"] = str(Path(config["OUTPUT_MODEL_PATH"]))
    if not config.get("OUTPUT_MODEL_NAME"):
        config["OUTPUT_MODEL_NAME"] = str(Path(config["OUTPUT_MODEL_PATH"]) / "model_classification.pt")
    if not config.get("ONNX_NAME"):
        config["ONNX_NAME"] = str(Path(config["OUTPUT_MODEL_PATH"]) / "model_classification.onnx")

    required_paths = ["TRAIN_DATASET", "VALIDATION_DATASET"]
    for key in required_paths:
        if not Path(config[key]).exists():
            raise FileNotFoundError(f"{key} file not found: {config[key]}")

    cast_config_values(config)
    return config


def cast_config_values(config):
    int_keys = ["RANDOM_SEED", "BATCH_SIZE", "EPOCHS", "PATIENCE", "HIDDEN_SIZE"]
    float_keys = ["LEARN_RATE", "WEIGHT_DEC"]
    bool_keys = ["TRIG_FUNCS", "USE_ONLY_CPU", "EXPORT_TO_ONNX"]

    for key in int_keys:
        config[key] = int(config[key])
    for key in float_keys:
        config[key] = float(config[key])
    for key in bool_keys:
        if isinstance(config[key], str):
            config[key] = str_to_bool(config[key])
        else:
            config[key] = bool(config[key])


if __name__ == "__main__":
    train_from_config(load_config(parse_args()))
