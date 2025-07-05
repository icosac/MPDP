import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.metrics import f1_score
from sklearn.ensemble import BaggingClassifier
import time


print("Reading ds")
clock = time.time()

data = pd.read_csv('/home/enrico/Projects/mpdp/small.csv', sep='\s+')
# data = pd.read_csv('/home/enrico/Projects/mpdp/big_smaller.csv', sep='\s+')
# data = pd.read_csv('/home/enrico/Projects/mpdp/big.csv', sep='\s+')

features = ['kmax', 'theta_i', 'theta_f', 'alpha_m', 'alpha_f']
target = 'id_man_comb'

X = data[features].values
y = data[target].values

print(f"Time to read ds: {time.time() - clock}")


print("Splitting ds")
clock = time.time()

reduce = False
if reduce:
    New_X, _, new_y, _ = train_test_split(X, y, test_size=0.999, shuffle=True)
else:
    New_X = X
    new_y = y

scaler = StandardScaler().fit(New_X)
New_X = scaler.fit_transform(New_X)

# Split the dataset into training (60%), validation (20%), and test (20%) sets
X_train, X_temp, y_train, y_temp = train_test_split(New_X, new_y, test_size=0.4, shuffle=True)
X_val, X_test, y_val, y_test = train_test_split(X_temp, y_temp, test_size=0.5, shuffle=True)

print(f"Time to split ds: {time.time() - clock}")
print(X_val.shape)


print("Training")
clock = time.time()
# Train SVM model
svm_model = SVC(kernel="rbf", probability=True)
# n_estimators = 2
# svm_model = BaggingClassifier(
#     SVC(
#         kernel='rbf', 
#         probability=True, 
#     ), 
#     max_samples=1.0 / n_estimators, 
#     n_estimators=n_estimators,
#     n_jobs=3
# )
svm_model.fit(X_train, y_train)
print(f"Time to train: {time.time() - clock}")


print("###############\nPredicting")
clock = time.time()
# Use the model to predict the validation set
# X_val = scaler.transform(X_val)
y_pred = svm_model.predict(X_val, )
print(f"Time to predict: {time.time() - clock}")
y_pred

# Calculate the accuracy of the model
accuracy = (y_pred == y_val).mean()
print(f"Accuracy: {accuracy:.2f}")
f_score = f1_score(y_val, y_pred, average='weighted')
print(f"F1 Score: {f_score:.2f}")


print("###############\nTesting")
clock = time.time()
# Use the model to predict the test set
# X_test = scaler.transform(X_test)
y_pred = svm_model.predict(X_test)
print(f"Time to predict: {time.time() - clock}")
y_pred

# Calculate the accuracy of the model
accuracy = (y_pred == y_test).mean()
print(f"Accuracy: {accuracy:.2f}")
f_score = f1_score(y_test, y_pred, average='weighted')
print(f"F1 Score: {f_score:.2f}")
