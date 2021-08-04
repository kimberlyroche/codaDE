import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel

# X must be samples x features
X = np.loadtxt("features.txt", delimiter=' ', skiprows=0, unpack=True)
X = X.transpose()
# y must be samples x 1 outcome
y = np.loadtxt("response_TPR.txt", delimiter=' ', skiprows=0, unpack=True)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

# RBF + noise + an offset
kernel = ConstantKernel(1.0) + ConstantKernel(1.0) * RBF(10) + WhiteKernel(5)
model = GaussianProcessRegressor(kernel=kernel)
model.fit(X_train, y_train)
# y_pred_tr, y_pred_tr_std = model.predict(X_train, return_std=True)
y_pred_te, y_pred_te_std = model.predict(X_test, return_std=True)

np.savetxt('pyGP_tests_TPR.txt', y_test, delimiter = ' ')
np.savetxt('pyGP_predictions_TPR.txt', y_pred_te, delimiter = ' ')
