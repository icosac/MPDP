import torch
from torch.utils.data import Dataset
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd

class DubinsDataset(Dataset):
    def __init__(self, data_path, transform=None):
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
        if all(col in self.data.columns for col in ['kmax', 'theta_i', 'theta_f', 'alpha_m', 'alpha_f', 'th_m']):
            # Extract features and target using column names
            self.features = self.data[['kmax', 'theta_i', 'theta_f', 'alpha_m', 'alpha_f']].values
            # Extract th_m and convert to sin and cos
            th_m_values = self.data['th_m'].values
            self.sin_th_m = np.sin(th_m_values)
            self.cos_th_m = np.cos(th_m_values)
            # extract the label of the maneuver
            self.labels = self.data['id_man_comb'].values
        else:
            print("Column names not found. Using positional columns instead.")
            # Extract features from the first 5 columns
            self.features = self.data.iloc[:, 0:5].values
            # Assuming th_m is in a specific column (adjust as needed)
            # Using column index 5 for th_m (adjust this based on your actual data structure)
            if self.data.shape[1] > 5:
                th_m_values = self.data.iloc[:, 5].values
                self.sin_th_m = np.sin(th_m_values)
                self.cos_th_m = np.cos(th_m_values)
                self.labels = self.data.iloc[:, 6].values
            else:
                raise ValueError("Cannot find 'th_m' column in the dataset")
        
        ### CLASSIFICATION TARGET ###
        
        # Normalize labels to start from 0
        unique_labels = np.unique(self.labels)
        print(f"Original unique labels: {unique_labels}")
        
        # Create a mapping from original labels to 0-indexed labels
        self.label_mapping = {original: idx for idx, original in enumerate(unique_labels)}
        self.inverse_mapping = {idx: original for original, idx in self.label_mapping.items()}
        
        # Map the labels to 0-indexed values
        self.labels = np.array([self.label_mapping[label] for label in self.labels])
        print(f"Remapped labels to range: 0-{len(self.label_mapping)-1}")
        print(f"Label mapping: {self.label_mapping}")
        
        ### REGRESSION TARGET ###
        
        # Combine sin and cos into targets array
        self.targets = np.column_stack((self.sin_th_m, self.cos_th_m))
        
        # Scale features
        self.scaler = StandardScaler()
        self.features = self.scaler.fit_transform(self.features)
        
        # No need to scale sin and cos as they're already normalized between -1 and 1
        
        # Convert to tensors
        self.features = torch.FloatTensor(self.features)
        self.targets = torch.FloatTensor(self.targets)
        self.labels = torch.LongTensor(self.labels)
        
        # Get number of unique classes
        self.num_classes = len(np.unique(self.labels))
        
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
            
        features = self.features[idx]
        target = self.targets[idx]
        label = self.labels[idx]
        
        if self.transform:
            features = self.transform(features)
            
        return features, target, label
    
    def get_scaler(self):
        return self.scaler
    
    def get_label_mapping(self):
        return self.label_mapping, self.inverse_mapping