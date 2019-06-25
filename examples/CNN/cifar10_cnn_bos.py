'''Train a simple deep CNN on the CIFAR10 small images dataset.

GPU run command with Theano backend (with TensorFlow, the GPU is automatically used):
    THEANO_FLAGS=mode=FAST_RUN,device=gpu,floatX=float32 python cifar10_cnn.py

It gets down to 0.65 test logloss in 25 epochs, and down to 0.55 after 50 epochs.
(it's still underfitting at that point, though).
'''

from __future__ import print_function
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Convolution2D, MaxPooling2D, BatchNormalization
from keras.utils import np_utils
from keras.optimizers import SGD
from keras.callbacks import ReduceLROnPlateau
from keras.backend.tensorflow_backend import set_session
from bos_function import run_BOS

import tensorflow as tf
import datetime as dt
import time
import numpy as np
import argparse
import scipy.io
import argparse
import keras

class BOS(keras.callbacks.Callback):
	
	def __init__(self, num_init = 20, incumbent = 0.5):
		self.val_epochs = []
		self.action_region = None
		self.grid_St = None
		self.num_init_curve = num_init
		self.incumbent = incumbent
		self.stop = 0
	
	def on_epoch_end(self, epoch, logs={}):
		scores = logs.get('val_acc')
		nb_epoch = self.params['nb_epoch']
		self.val_epochs.append(scores)
		if (epoch+1 == self.num_init_curve):
			self.action_region, self.grid_St = run_BOS(1 - np.array(self.val_epochs), 
					self.incumbent, self.params['nb_epoch'])

		if (epoch >= self.num_init_curve) and (epoch < nb_epoch - 1):
				state = np.sum(1 - np.array(self.val_epochs[self.num_init_curve:])) / (epoch - self.num_init_curve + 1)
				ind_state = np.max(np.nonzero(state > self.grid_St)[0])
				
				action_to_take = self.action_region[epoch - self.num_init_curve, ind_state]
				if (action_to_take) == 1 or (scores >= self.incumbent):
					self.stop = 1
					self.model.stop_training = True
					
				elif action_to_take == 2:
					self.stop = -1
					self.model.stop_training = True
		
		elif (epoch == nb_epoch-1):
				self.stop = -1

def load_data(method='full'):

	if method == 'full':
		x_train = np.load('./cifar10/xtrain.npy')
		y_train = np.load('./cifar10/ytrain.npy')
	else:
		x_train = np.load('./cifar10/xtrain_small.npy')
		y_train = np.load('./cifar10/ytrain_small.npy')

	x_test = np.load('./cifar10/xval.npy')
	y_test = np.load('./cifar10/yval.npy')
    
	return x_train, y_train, x_test, y_test

def cnn_cifar10_train(X_train, y_train, X_test, y_test, epochs, dropout, lr, batch_size, method):
	
	config = tf.ConfigProto()
	config.gpu_options.allow_growth = True
	set_session(tf.Session(config=config))
	
	start_time = dt.datetime.now()

	#batch_size = b_size
	nb_classes = 10
	nb_epoch = epochs
	data_augmentation = True
	num_init_curve = 20
	incumbent = 0.5

	# input image dimensions
	img_rows, img_cols = 32, 32
	# The CIFAR10 images are RGB.
	img_channels = 3

	# Convert class vectors to binary class matrices.
	Y_train = np_utils.to_categorical(y_train, nb_classes)
	Y_test = np_utils.to_categorical(y_test, nb_classes)

	model = Sequential()

	model.add(Convolution2D(32, 3, 3, border_mode='same',
							input_shape=X_train.shape[1:]))
	model.add(Activation('relu'))
	model.add(Convolution2D(32, 3, 3))
	model.add(Activation('relu'))
	model.add(MaxPooling2D(pool_size=(2, 2)))
	model.add(Dropout(dropout[0]))

	model.add(Convolution2D(64, 3, 3, border_mode='same'))
	model.add(Activation('relu'))
	model.add(Convolution2D(64, 3, 3))
	model.add(Activation('relu'))
	model.add(MaxPooling2D(pool_size=(2, 2)))
	model.add(Dropout(dropout[1]))

	model.add(Flatten())
	model.add(Dense(512))
	model.add(Activation('relu'))
	model.add(Dropout(dropout[2]))

	model.add(Dense(nb_classes))
	model.add(Activation('softmax'))

	# Let's train the model using RMSprop
	sgd = SGD(lr=lr, decay=1e-5, momentum=0.9, nesterov=True)
	model.compile(loss='categorical_crossentropy',
				  optimizer=sgd,
				  metrics=['accuracy'])

	X_train = X_train.astype('float32')
	X_test = X_test.astype('float32')
	X_train /= 255
	X_test /= 255
	
	# Create callbacks
	lr_change = ReduceLROnPlateau(monitor='val_loss', factor=0.2,
					  patience=50, min_lr=0.000001, verbose=1)
	call_bos = BOS(num_init=num_init_curve, incumbent=incumbent)

	if not data_augmentation:
		print('Not using data augmentation.')
		model.fit(X_train, Y_train,
				  batch_size=batch_size,
				  nb_epoch=nb_epoch,
				  validation_data=(X_test, Y_test),
				  shuffle=True, verbose=2)
	else:
		print('Using real-time data augmentation.')
		# This will do preprocessing and realtime data augmentation:
		datagen = ImageDataGenerator(
			featurewise_center=False,  # set input mean to 0 over the dataset
			samplewise_center=False,  # set each sample mean to 0
			featurewise_std_normalization=False,  # divide inputs by std of the dataset
			samplewise_std_normalization=False,  # divide each input by its std
			zca_whitening=False,  # apply ZCA whitening
			rotation_range=0,  # randomly rotate images in the range (degrees, 0 to 180)
			width_shift_range=0.1,  # randomly shift images horizontally (fraction of total width)
			height_shift_range=0.1,  # randomly shift images vertically (fraction of total height)
			horizontal_flip=True,  # randomly flip images
			vertical_flip=False)  # randomly flip images

		# Compute quantities required for featurewise normalization
		# (std, mean, and principal components if ZCA whitening is applied).
		datagen.fit(X_train)

		# Fit the model on the batches generated by datagen.flow().
		if method == 'full':
			hist = model.fit_generator(datagen.flow(X_train, Y_train,
							batch_size=batch_size),
						samples_per_epoch=X_train.shape[0],
						nb_epoch=nb_epoch,
						validation_data=(X_test, Y_test), verbose=2,
						callbacks=[lr_change])
			acc = hist.history['val_acc'][nb_epoch-1]
		else:
			hist = model.fit_generator(datagen.flow(X_train, Y_train,
							batch_size=batch_size),
						samples_per_epoch=X_train.shape[0],
						nb_epoch=nb_epoch,
						validation_data=(X_test, Y_test), verbose=2,
						callbacks=[lr_change, call_bos])
			acc = call_bos.stop
			
	end_time = dt.datetime.now()
	elapsed_time= (end_time - start_time).seconds/3600.

	return acc, elapsed_time

if __name__=="__main__":
	
	newParser = argparse.ArgumentParser()
	
	newParser.add_argument("--method", dest="TRAIN_METHOD", default='full')
	newParser.add_argument("--lr", dest="lr", default=0.01, type=float)
	newParser.add_argument("--d1", dest="d1", default=0.5, type=float)
	newParser.add_argument("--d2", dest="d2", default=0.5, type=float)
	newParser.add_argument("--d3", dest="d3", default=0.5, type=float)
	newParser.add_argument("--batch_size", dest="BATCH_SIZE", default=128, type=int)
	newParser.add_argument("--epochs", dest="EPOCHS", default=1000, type=int)
	newParser.add_argument("--output", dest="fout", default="score_time.txt")
	args = newParser.parse_args()
	
	epochs = args.EPOCHS
	dropout = [args.d1, args.d2, args.d3]
	batch_size = args.BATCH_SIZE
	lr = args.lr
	file_out = args.fout

	print("loading data...")
	x_train, y_train, x_test, y_test = load_data(method=args.TRAIN_METHOD)

	print(epochs, dropout, lr, batch_size)
	score, elapsed_time = cnn_cifar10_train(x_train, y_train, x_test, y_test,
			epochs, dropout, lr, batch_size, args.TRAIN_METHOD)

	fpt = open(file_out, 'w')
	fpt.write(str(score)+" "+str(elapsed_time))
	fpt.close()
	print('epochs: ', epochs, ' dropout: ', dropout, ' lr: ', lr, ' batch_size: ', batch_size)
	print('performance: ', score, 'time: ', elapsed_time, '\n')
