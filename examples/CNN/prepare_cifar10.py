# Prepare the dataet for CIFAR10_CNN_BOS
from keras.datasets import cifar10
import numpy as np
import os

if not os.path.isdir("./cifar10"):
	os.mkdir("./cifar10")


(x_train, y_train), (x_val, y_val) = cifar10.load_data()

np.save('./cifar10/yval', y_val)
np.save('./cifar10/xval', x_val)

np.save('./cifar10/xtrain', x_train)
np.save('./cifar10/ytrain', y_train)

num = 10000	# no. of training data for the auxiliary function
ind = np.arange(50000)

np.random.shuffle(ind)
xs = x_train[ind[0:num]]
ys = y_train[ind[0:num]]

np.save('./cifar10/xtrain_small', xs)
np.save('./cifar10/ytrain_small', ys)

