import torch
from torch.utils.data import Dataset
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd

class DubinsDatasetRectangle(Dataset):
    def __init__(self, data_path, transform=None, use_trigonometric_features=False):
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
        if all(col in self.data.columns for col in ['kmax', 'xi', 'xm', 'ym', 'xf', 'theta_i', 'theta_f', 'th_m']):
            if use_trigonometric_features:
                print("Using trigonometric features for 'theta_i' and 'theta_f'.")
                # Extract features and target using column names
                th_is = self.data['theta_i'].values
                th_fs = self.data['theta_f'].values
                
                sin_th_is = np.sin(th_is)
                cos_th_is = np.cos(th_is)
                sin_th_fs = np.sin(th_fs)
                cos_th_fs = np.cos(th_fs)

                self.features = np.column_stack((self.data['kmax'].values, self.data['xi'].values, 
                                                 self.data['xm'].values, self.data['ym'].values, 
                                                 self.data['xf'].values, sin_th_is, cos_th_is, 
                                                 sin_th_fs, cos_th_fs))
                # self.features = self.data[['kmax', 'theta_i', 'theta_f', 'alpha_m', 'alpha_f']].values
                # Extract th_m and convert to sin and cos
                th_m_values = self.data['th_m'].values
                self.sin_th_m = np.sin(th_m_values)
                self.cos_th_m = np.cos(th_m_values)
            else:
                # Extract features using column names without trigonometric transformations
                self.features = self.data[['kmax', 'xi', 'xm', 'ym', 'xf', 'theta_i', 'theta_f']].values

            # Extract th_m and convert to sin and cos
            th_m_values = self.data['th_m'].values
            self.sin_th_m = np.sin(th_m_values)
            self.cos_th_m = np.cos(th_m_values)
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
            else:
                raise ValueError("Cannot find 'th_m' column in the dataset")
        
        # Combine sin and cos into targets array
        self.targets = np.column_stack((self.sin_th_m, self.cos_th_m))
        
        # Scale features
        self.scaler = StandardScaler()
        self.features = self.scaler.fit_transform(self.features)
        
        # No need to scale sin and cos as they're already normalized between -1 and 1
        
        # Convert to tensors
        self.features = torch.FloatTensor(self.features)
        self.targets = torch.FloatTensor(self.targets)
        
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
            
        features = self.features[idx]
        target = self.targets[idx]
        
        if self.transform:
            features = self.transform(features)
            
        return features, target
    
    def get_scaler(self):
        return self.scaler
    
    def get_num_features(self):
        return self.features.shape[1]