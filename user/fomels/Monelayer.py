#!/usr/bin/env python
""" One layer NN training.

Reproducing the training process of the following example:
https://github.com/seg/tutorials-2018/tree/master/1808_Neural_networks
This is a minimal one-layer network:
    x -> a -> z == y
The training uses L2 loss function.
"""
from __future__ import print_function
import sys
try:
    import numpy as np
    from numpy.random import randn
    import rsf.api as rsf
except Exception as e:
    print('ERROR: need numpy')
    sys.exit(1)

# Network helper functions.
def sigma(z, forward=True):
    if forward:
        return 1 / (1 + np.exp(-z))
    else:
        return z * (1 - z)
def tanh(z, forward=True):
    if forward:
        return (np.exp(z) - np.exp(-z)) / (np.exp(z) + np.exp(-z)) 
    else:
        return 1 - tanh(z)**2
def relu(z, forward=True):
    if forward:
        return z * (z > 0) 
    else:
        return 1 * (z > 0)

def forward(xi, W1, b1, W2, b2):
    z1 = W1.dot(xi) + b1
    a1 = sigma(z1)
    z2 = W2.dot(a1) + b2
    return z2, a1

def backward(xi, yi,
             a1, z2,
             params,
             learning_rate):

    err_output = z2 - yi
    grad_W2 = err_output * a1
    params['W2'] -= learning_rate * grad_W2

    grad_b2 = err_output
    params['b2'] -= learning_rate * grad_b2

    derivative = sigma(a1, forward=False)
    err_hidden = err_output * derivative * params['W2']
    grad_W1 = err_hidden[:, None].dot(xi[None, :])
    params['W1'] -= learning_rate * grad_W1
    
    grad_b1 = err_hidden
    params['b1'] -= learning_rate * grad_b1
    
    return params

def initialize_params(units, features):
    np.random.seed(42)
    params = {
        "W1": 0.1 * randn(units, features),
        "b1": np.zeros(shape=units),

        "W2": 0.1 * randn(units),
        "b2": np.zeros(shape=1)
    }
    return params


# Get input parameters and 'source' files. These are RSF objects that works as
# pointers to the actual data/file.
par_obj = rsf.Par()
x_train_obj = rsf.Input() # first input as default
y_train_obj = rsf.Input('ytrain')
x_val_obj = rsf.Input('xval')
y_val_obj = rsf.Input('yval')

# Assign output 'target' files.
loss_train_obj = rsf.Output() # first output as default
loss_val_obj = rsf.Output('lossval')

# Parse parameter values from input object.
num_epochs = par_obj.int('nepochs')
learning_rate = par_obj.float('lr')
hidden_units = par_obj.int('hidden') # number of hidden units

# Get size of the input arrays and load them. (Initialize + Read)
# Initialize numpy arrays for training data.
n1 = x_train_obj.int('n1') # for x_train, n1==data sample length
n2 = x_train_obj.int('n2') # for x_train, n2==number of features
x_train = np.zeros((n2,n1), np.single) # shape==(nfeature x ndata)
y_train = np.zeros((1,n1), np.single) # shape==(1 x ndata)
# Similar for val data.
n1 = x_val_obj.int('n1')
n2 = x_val_obj.int('n2')
x_val = np.zeros((n2,n1), np.single)
y_val = np.zeros((1,n1), np.single)
# Read from objects. If dimension mismatches, prog will exit.
x_train_obj.read(x_train)
y_train_obj.read(y_train)
x_val_obj.read(x_val)
y_val_obj.read(y_val)


# Arrange the dimensions of the input arrays.
x_train, x_val = x_train.T, x_val.T # now (ndata x nfeature)
assert x_train.shape[-1] == x_val.shape[-1] # number of features should be equal
y_train = y_train.reshape(-1) # squeeze from (1 x ndata) to (ndata)
y_val = y_val.reshape(-1)

# Training process as below, similar to the original post.
net_params = initialize_params(hidden_units, x_train.shape[-1])
data_train = list(zip(x_train, y_train))
data_val = list(zip(x_val, y_val))
loss_train_history, loss_val_history = [], []
for i in range(num_epochs):
    # Validation. We do this first to get the same training state as for the
    # training data (below).
    np.random.shuffle(data_val)
    loss = 0.0
    for x, y in data_val:
        z, a = forward(x, **net_params)
        loss += np.square(z - y)
    loss_val_history.append(loss / y_val.size)

    # Training.
    np.random.shuffle(data_train)
    loss = 0.0
    for x, y in data_train:
        z, a = forward(x, **net_params)
        net_params = backward(x, y, a, z, net_params, learning_rate)
        loss += np.square(z - y)
    loss_train_history.append(loss / y_train.size)


# Output history logs to files.
loss_train = np.array(loss_train_history)
loss_val = np.array(loss_val_history)
loss_train_obj.put('n1', num_epochs)
loss_train_obj.put('o1', 1)
loss_train_obj.put('d1', 1)
loss_train_obj.put('label1', "Epoch")
loss_train_obj.put('unit1', "")
loss_train_obj.put('n2', 1)
loss_train_obj.put('label2', "Loss")
loss_train_obj.put('unit2', "")
loss_val_obj.put('n1', num_epochs)
loss_val_obj.put('o1', 1)
loss_val_obj.put('d1', 1)
loss_val_obj.put('label1', "Epoch")
loss_val_obj.put('unit1', "")
loss_val_obj.put('n2', 1)
loss_val_obj.put('label2', "Loss")
loss_val_obj.put('unit2', "")
loss_train_obj.write(loss_train)
loss_val_obj.write(loss_val)
loss_train_obj.close()
loss_val_obj.close()
