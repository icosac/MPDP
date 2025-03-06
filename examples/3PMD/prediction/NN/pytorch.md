# Correct installation

Follow this [link](https://pytorch.org/get-started/locally/) as the correct packages may depend on 
your installation. 

## Linux

### No CUDA

```shell
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
```

### CUDA 11.8

```shell
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
```

### CUDA 12.4

```shell
pip3 install torch torchvision torchaudio
```

### CUDA 12.6

```shell
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu126
```

## MacOS

```shell
pip3 install torch torchvision torchaudio
```