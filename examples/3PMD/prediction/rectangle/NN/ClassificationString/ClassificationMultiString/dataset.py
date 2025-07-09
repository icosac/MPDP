import torch
from torch.utils.data import Dataset
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd

class DubinsDatasetRectangle(Dataset):
    def __init__(self, data_path, transform=None, use_trigonometric_features=True):
        """
        Args:
            data_path (string): Path to the CSV file with annotations.
            transform (callable, optional): Optional transform to be applied on a sample.
        """
        # Load data from CSV
        self.data = pd.read_csv(data_path, delim_whitespace=True)  # Use whitespace as delimiter
        print(f"Loaded columns: {self.data.columns.tolist()}")
        self.transform = transform
        
        # Check if columns exist, otherwise use positional columns
        if all(col in self.data.columns for col in ['kmax', 'xi', 'xm', 'ym', 'xf', 'theta_i', 'theta_f']):
            # Extract features and labels using column names

            if use_trigonometric_features:
                print("Using trigonometric features.")
                theta_is = self.data['theta_i'].values
                theta_fs = self.data['theta_f'].values
                
                sin_theta_is = np.sin(theta_is)
                cos_theta_is = np.cos(theta_is)
                sin_theta_fs = np.sin(theta_fs)
                cos_theta_fs = np.cos(theta_fs)

                self.features = np.column_stack((self.data['kmax'].values, self.data['xi'].values, self.data['xm'].values,
                                                self.data['ym'].values, self.data['xf'].values, sin_theta_is,
                                                cos_theta_is, sin_theta_fs, cos_theta_fs))
            else:
                self.features = self.data[['kmax', 'xi', 'xm', 'ym', 'xf', 'theta_i', 'theta_f']].values

            self.labels = self.data['id_man_comb'].values
        else:
            print("Column names not found. Using positional columns instead.")
            # Extract features from the first 7 columns and labels from the 9th column
            self.features = self.data.iloc[:, 0:7].values
            self.labels = self.data.iloc[:, 8].values  
        
        # Convert all labels to string type
        self.labels = self.labels.astype(str)
        unique_labels = np.unique(self.labels)
        print(f"Original unique labels: {unique_labels}")
        # Create a mapping from string labels to 0-indexed integers
        self.label_mapping = {original: idx for idx, original in enumerate(unique_labels)}
        self.inverse_mapping = {idx: original for original, idx in self.label_mapping.items()}
        # Map the string labels to 0-indexed integer values for training
        self.labels_idx = np.array([self.label_mapping[label] for label in self.labels])
        print(f"Label mapping: {self.label_mapping}")
        # Scale features
        self.scaler = StandardScaler()
        self.features = self.scaler.fit_transform(self.features)
        # Convert to tensors
        self.features = torch.FloatTensor(self.features)
        self.labels_idx = torch.LongTensor(self.labels_idx)
        # Get number of unique classes
        self.num_classes = len(unique_labels)
        
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
        features = self.features[idx]
        label_idx = self.labels_idx[idx]
        if self.transform:
            features = self.transform(features)
        return features, label_idx
    
    def get_scaler(self):
        return self.scaler
        
    def get_label_mapping(self):
        return self.label_mapping, self.inverse_mapping

    def get_num_features(self):
        return self.features.shape[1]
