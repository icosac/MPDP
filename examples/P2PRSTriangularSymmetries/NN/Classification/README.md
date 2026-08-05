# RS Triangular Maneuver Classification

This example trains a small PyTorch classifier for the point-to-point Reeds-Shepp
triangular-symmetry dataset.

The expected dataset columns are:

```text
thi thf kmax len man man_str n_seg
```

Only `thi`, `thf`, and `kmax` are used as inputs. The label is `man`. With
`TRIG_FUNCS: true`, the network is trained on:

```text
sin(thi), cos(thi), sin(thf), cos(thf), kmax
```

## Train

Place the datasets here:

```text
examples/P2PRSTriangularSymmetries/datasets/DS_RS_triangular_training.csv
examples/P2PRSTriangularSymmetries/datasets/DS_RS_triangular_validation.csv
```

Then run:

```bash
cd examples/P2PRSTriangularSymmetries/NN/Classification
python3 main.py --yaml-config config/rs_triangular_class_config.yaml
```

You can also override paths directly:

```bash
python3 main.py \
  --train-dataset /path/to/DS_RS_triangular_training.csv \
  --validation-dataset /path/to/DS_RS_triangular_validation.csv
```

Training writes the checkpoint, scaler, optional ONNX model, and plots under
`models/` and `plots/`.

## Inference

For a single raw input in `thi,thf,kmax` order:

```bash
python3 inference.py --model models/model_classification.pt --input "0.0087509544668126,0.0087509544668126,0.1"
```

or with the YAML config:

```bash
python3 inference.py --yaml-config rs_triangular_inf_config.yaml
```

For a CSV with the same columns as the training data:

```bash
python3 inference.py --model models/model_classification.pt --input /path/to/DS_RS_triangular_validation.csv
```

## C++ ONNX Inference

CPU build:

```bash
cd examples/P2PRSTriangularSymmetries/NN/Classification/inference_cpp
cmake -S . -B build
cmake --build build -j
```

CUDA build on Linux:

```bash
cmake -S . -B build-cuda -DRS_TRIANGULAR_USE_CUDA=ON
cmake --build build-cuda -j
```

Run CUDA inference with:

```bash
./build-cuda/rs_triangular_inference ../models/model_classification.onnx "0.0087509544668126,0.0087509544668126,0.1" --cuda
```

Use `--cuda-device N` to select a GPU. If you pass `ONNXRUNTIME_ROOT`
manually, it must point to an extracted ONNX Runtime GPU package containing
`lib/libonnxruntime_providers_cuda.so`.
