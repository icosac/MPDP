import numpy as np
import pandas as pd
import torch
from sklearn.preprocessing import StandardScaler
from torch.utils.data import Dataset


RAW_FEATURE_COLUMNS = ["thi", "thf", "kmax"]
LABEL_COLUMN = "man"


class RSTriangularDataset(Dataset):
    """Dataset for point-to-point Reeds-Shepp triangular-symmetry classification."""

    def __init__(
        self,
        data_path,
        scaler=None,
        label_mapping=None,
        use_trigonometric_features=True,
        fit_scaler=False,
    ):
        self.data_path = data_path
        self.data = pd.read_csv(data_path, sep=r"[\s,]+", engine="python")
        self.use_trigonometric_features = use_trigonometric_features

        missing = [col for col in RAW_FEATURE_COLUMNS + [LABEL_COLUMN] if col not in self.data.columns]
        if missing:
            raise ValueError(
                f"Dataset {data_path} is missing required columns: {missing}. "
                f"Available columns: {self.data.columns.tolist()}"
            )

        raw_features = self.data[RAW_FEATURE_COLUMNS].to_numpy(dtype=np.float32)
        features = build_features(raw_features, use_trigonometric_features)
        labels = self.data[LABEL_COLUMN].to_numpy(dtype=np.int64)

        if label_mapping is None:
            unique_labels = sorted(np.unique(labels).tolist())
            self.label_mapping = {int(original): idx for idx, original in enumerate(unique_labels)}
        else:
            self.label_mapping = {int(original): int(idx) for original, idx in label_mapping.items()}

        unknown_labels = sorted(set(labels.tolist()) - set(self.label_mapping.keys()))
        if unknown_labels:
            raise ValueError(
                f"Dataset {data_path} contains labels not present in the training mapping: {unknown_labels}"
            )

        self.inverse_mapping = {idx: original for original, idx in self.label_mapping.items()}
        mapped_labels = np.array([self.label_mapping[int(label)] for label in labels], dtype=np.int64)

        if scaler is None:
            self.scaler = StandardScaler()
            if fit_scaler:
                features = self.scaler.fit_transform(features)
            else:
                features = self.scaler.fit_transform(features)
        else:
            self.scaler = scaler
            features = self.scaler.transform(features)

        self.raw_features = torch.as_tensor(raw_features, dtype=torch.float32)
        self.features = torch.as_tensor(features, dtype=torch.float32)
        self.labels = torch.as_tensor(mapped_labels, dtype=torch.long)
        self.num_classes = len(self.label_mapping)

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        return self.features[idx], self.labels[idx]

    def get_scaler(self):
        return self.scaler

    def get_label_mapping(self):
        return self.label_mapping, self.inverse_mapping

    def get_num_features(self):
        return self.features.shape[1]

    def get_feature_names(self):
        if self.use_trigonometric_features:
            return ["sin_thi", "cos_thi", "sin_thf", "cos_thf", "kmax"]
        return RAW_FEATURE_COLUMNS.copy()


def build_features(raw_features, use_trigonometric_features=True):
    """Build model features from raw [thi, thf, kmax] rows."""
    raw_features = np.asarray(raw_features, dtype=np.float32)
    if raw_features.ndim == 1:
        raw_features = raw_features.reshape(1, -1)
    if raw_features.shape[1] != len(RAW_FEATURE_COLUMNS):
        raise ValueError(
            f"Expected raw features with columns {RAW_FEATURE_COLUMNS}, got shape {raw_features.shape}"
        )

    thi = raw_features[:, 0]
    thf = raw_features[:, 1]
    kmax = raw_features[:, 2]

    if not use_trigonometric_features:
        return raw_features.astype(np.float32)

    return np.column_stack(
        (
            np.sin(thi),
            np.cos(thi),
            np.sin(thf),
            np.cos(thf),
            kmax,
        )
    ).astype(np.float32)
