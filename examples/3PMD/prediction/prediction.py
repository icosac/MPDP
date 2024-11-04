import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.optimizers import Adam

# Load the dataset
data = pd.read_csv('/Users/enrico/Projects/mpdp/3PDS_Circle32_1_1_1.csv', sep='\s+')

# Separate features and target variables
features = ['kmax', 'theta_i', 'theta_f', 'alpha_m', 'alpha_f', 'len']
targets = ['id_man_comb', 'th_m']

X = data[features].values
y = data[targets].values

# Split the dataset into training (60%), validation (20%), and test (20%) sets
X_train, X_temp, y_train, y_temp = train_test_split(X, y, test_size=0.4, random_state=42)
X_val, X_test, y_val, y_test = train_test_split(X_temp, y_temp, test_size=0.5, random_state=42)

# Standardize the features
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_val = scaler.transform(X_val)
X_test = scaler.transform(X_test)

# Build the neural network model
model = Sequential([
    Dense(64, activation='relu', input_shape=(X_train.shape[1],)),
    Dense(32, activation='relu'),
    Dense(16, activation='relu'),
    Dense(2)  # Output layer with 2 outputs for multi-output regression
])

# Compile the model
model.compile(optimizer=Adam(learning_rate=0.001), loss='mse')

# Train the model
history = model.fit(X_train, y_train, validation_data=(X_val, y_val), epochs=10, batch_size=32, verbose=1)

# Evaluate the model on the test set
test_loss = model.evaluate(X_test, y_test)
print(f'Test Mean Squared Error: {test_loss}')

# Plot training and validation loss over epochs
plt.figure(figsize=(10, 6))
plt.plot(history.history['loss'], label='Training Loss')
plt.plot(history.history['val_loss'], label='Validation Loss')
plt.title('Training and Validation Loss')
plt.xlabel('Epochs')
plt.ylabel('Mean Squared Error')
plt.legend()
plt.show()

# Function to make predictions
def predict_outputs(kmax, theta_i, theta_f, alpha_m, alpha_f, length):
    input_data = scaler.transform([[kmax, theta_i, theta_f, alpha_m, alpha_f, length]])
    prediction = model.predict(input_data)
    id_pred = prediction[0][0]  # Predicted id_man_comb
    double_pred = prediction[0][1]  # Predicted target_double
    return id_pred, double_pred

# Example prediction
example_prediction_id, example_prediction_double = predict_outputs(1.5, 0.3, 0.7, 0.2, 0.4, 5.0)
print(f'Predicted id_man_comb: {example_prediction_id}')
print(f'Predicted target_double: {example_prediction_double}')
