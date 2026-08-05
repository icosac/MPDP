import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import torch
import yaml

from dataset import LABEL_COLUMN, RAW_FEATURE_COLUMNS, build_features
from model import NeuralNet


def load_inference_model(model_path, device):
    checkpoint = torch.load(model_path, map_location=device)
    model = NeuralNet(
        checkpoint["input_size"],
        checkpoint["hidden_size"],
        checkpoint["num_classes"],
        layers=checkpoint["layers"],
    )
    model.load_state_dict(checkpoint["model_state_dict"])
    model.to(device)
    model.eval()

    scaler_path = Path(model_path).with_suffix(".scaler.joblib")
    scaler = joblib.load(scaler_path)
    return model, scaler, checkpoint


def predict(model, scaler, checkpoint, raw_features, device, topk=3):
    features = build_features(raw_features, checkpoint["use_trigonometric_features"])
    features = scaler.transform(features)
    tensor = torch.as_tensor(features, dtype=torch.float32, device=device)

    with torch.no_grad():
        logits = model(tensor)
        probs = torch.softmax(logits, dim=1)
        top_probs, top_indices = torch.topk(probs, min(topk, checkpoint["num_classes"]), dim=1)

    inverse_mapping = {int(k): int(v) for k, v in checkpoint["inverse_mapping"].items()}
    top_indices_np = top_indices.cpu().numpy()
    top_mans = np.vectorize(inverse_mapping.get)(top_indices_np)
    return top_mans, top_probs.cpu().numpy()


def main():
    parser = argparse.ArgumentParser(description="Run inference with a trained RS triangular classifier.")
    parser.add_argument("--yaml-config", type=str, default=None, help="Path to inference YAML configuration.")
    parser.add_argument("--model", type=str, default="models/model_classification.pt")
    parser.add_argument("--input", type=str, default=None, help="Either 'thi,thf,kmax' or a whitespace CSV path.")
    parser.add_argument("--topk", type=int, default=3)
    parser.add_argument("--use-only-cpu", action="store_true")
    args = parser.parse_args()

    config = {}
    if args.yaml_config:
        with open(args.yaml_config, "r", encoding="utf-8") as stream:
            config = yaml.safe_load(stream) or {}

    model_path = args.model if args.model != "models/model_classification.pt" else config.get("MODEL", args.model)
    topk = args.topk if args.topk != 3 else int(config.get("TOPK", args.topk))
    use_only_cpu = args.use_only_cpu or bool(config.get("USE_ONLY_CPU", False))
    raw_input = args.input or config.get("INPUT")
    if raw_input is None:
        raise ValueError("Provide --input or INPUT in the YAML config.")

    device = torch.device("cuda:0" if (not use_only_cpu and torch.cuda.is_available()) else "cpu")
    model, scaler, checkpoint = load_inference_model(model_path, device)

    input_path = Path(raw_input)
    if input_path.exists():
        data = pd.read_csv(input_path, sep=r"[\s,]+", engine="python")
        raw_features = data[RAW_FEATURE_COLUMNS].to_numpy(dtype=np.float32)
        top_mans, top_probs = predict(model, scaler, checkpoint, raw_features, device, topk)

        correct = 0
        has_labels = LABEL_COLUMN in data.columns
        for idx, (mans, probs) in enumerate(zip(top_mans, top_probs)):
            print(f"Sample {idx}: top-{len(mans)} man {mans.tolist()} probabilities {probs.tolist()}")
            if has_labels and int(data.iloc[idx][LABEL_COLUMN]) in mans:
                correct += 1

        if has_labels:
            print(f"Top-{topk} accuracy: {correct / len(data):.6f}")
    else:
        raw_features = np.array([float(value) for value in raw_input.split(",")], dtype=np.float32)
        top_mans, top_probs = predict(model, scaler, checkpoint, raw_features, device, topk)
        print(f"Top-{len(top_mans[0])} man: {top_mans[0].tolist()}")
        print(f"Probabilities: {top_probs[0].tolist()}")


if __name__ == "__main__":
    main()
