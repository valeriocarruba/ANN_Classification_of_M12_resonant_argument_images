# -*- coding: utf-8 -*-
"""
# ===========================================================================
# ===========================================================================
# !==   Safwan ALJBAAE, Valerio Carruba                                    ==
# !==   November 2020                                                      ==
# ===========================================================================
"""

import tensorflow as tf
import matplotlib.pyplot as plt
from tensorflow import keras
from keras.utils.vis_utils import plot_model
from keras.callbacks import ModelCheckpoint
import numpy as np
import pandas as pd
from PIL import Image
from sys import getsizeof
import copy

class_names = ['circulation state', 'switching orbits', 'libration state']

data_ast = pd.read_csv('../../ALL_PNG_varpic/all_m12_status_new',
                       skiprows=0,
                       header=None,
                       delim_whitespace=True,
                       index_col=None,
                       names=['id', 'a', 'e', 'sin_i', 'Mag', 'label'],
                       low_memory=False,
                       dtype={'id': np.integer,
                              'a': np.float64,
                              'e': np.float64,
                              'sin_i': np.float64,
                              'Mag': np.float64,
                              'label': np.integer
                              }
                       )
n_train = int(len(data_ast))
train_data = data_ast.iloc[0:n_train, :]
train_id = train_data.id
img = ['../../ALL_PNG_varpic/fig_res_' + str("{:07d}".format(ast_id)) + '.png' for ast_id in train_id]
width, height = Image.open(img[0]).convert('1').size
train_images = [np.array(Image.open(x).convert('1').getdata()).reshape(width, height) for x in img]
train_label = train_data.label

print(f'The size of the variable train_images is : {getsizeof(train_images)} bytes')
print(f'We have {len(train_images)} images ({width} X {height} pixels) in the training set, belonging to '
      f'{len(set(train_label))} classes:')
for i in range(len(set(train_label))):
    print(f'   {len([x for x in train_label if x == i])} asteroids in {class_names[i]} (label: {i})')
print()

test_data = pd.read_csv('m12_status',
                        skiprows=0,
                        header=None,
                        delim_whitespace=True,
                        index_col=None,
                        names=['id', 'a', 'e', 'sin_i', 'Mag', 'label'],
                        low_memory=False,
                        dtype={'id': np.integer,
                               'a': np.float64,
                               'e': np.float64,
                               'sin_i': np.float64,
                               'Mag': np.float64}
                        )
# n_test=int(len(test_data))
# test_data= test_data.iloc[0:n_test, :]
test_id = list(test_data.id)
img = ['./fig_res_' + str("{:07d}".format(ast_id)) + '.png' for ast_id in test_id]
width, height = Image.open(img[0]).convert('1').size
test_images = [np.array(Image.open(x).convert('1').getdata()).reshape(width, height) for x in img]

print(f'The size of the variable test_images is: {getsizeof(train_images)} bytes')

min_pixel = min(list(map(min, train_images[0])))
max_pixel = max(list(map(max, train_images[0])))
print(f'The pixel values of each image varies from {min_pixel} to {max_pixel}')

# preprocessing the data: rescale the pixels value to range from 0 to 1
train_images = train_images / max_pixel
test_images = test_images / max_pixel

# Set up the model
model = keras.Sequential([
    keras.layers.Flatten(input_shape=(train_images.shape[1], train_images.shape[2])),
    keras.layers.Dense(128, activation=tf.nn.relu),
    keras.layers.Dense(64, activation=tf.nn.relu),
    keras.layers.Dense(len(class_names), activation=tf.nn.softmax)
])
model.summary()
plot_model(model, to_file='model_plot.png', show_shapes=True, show_layer_names=True)
# compile the model
model.compile(optimizer='adam',
              loss='sparse_categorical_crossentropy',
              metrics='accuracy')

# This checkpoint object will store the model parameters
# in the file "weights.hdf5"
checkpoint = ModelCheckpoint('./weights.hdf5',
                             save_weights_only=True,
                             monitor='accuracy',
                             mode='max',
                             verbose=1,
                             save_best_only=True,)

# fit the model, using the checkpoint as a callback
x = model.fit(train_images, train_label, epochs=30, callbacks=[checkpoint], verbose=0)
model.load_weights('./weights.hdf5')

# Plot the history of the model
fig = plt.figure()
figure = fig.add_subplot(111)
figure.plot(x.epoch, x.history['accuracy'])
plt.xlabel('Epoch')
plt.ylabel('Accuracy')

fig.savefig('history_model.png', format='png', dpi=300)
plt.close(fig)

predictions = model.predict(test_images)
predict_label = [int(np.argmax(x)) for x in predictions]
predict_acc = [100 * max(x) for x in predictions]

predicted_data = copy.deepcopy(test_data)
predicted_data['predicted_label'] = list(predict_label)
predicted_data.to_csv(r'predicted_data.csv', index=False, header=False, sep=' ', float_format='%.7f')
print()

# Compute the metrics values for the predicted test set
#for type_orbit in set(test_label):
#    data = predicted_data[predicted_data.predicted_label == type_orbit]
#    TP = len(data[data.label == type_orbit])
#    FP = len(data) - TP
#    data = predicted_data[predicted_data.label == type_orbit]
#    TN = len(data[data.predicted_label == type_orbit])
#    FN = len(data) - TP
#    purity = float(TP)/float(TP+FP)
#    compl = float(TP)/float(TP+FN)
#    merit = np.sqrt(4*purity**2+compl**2)/np.sqrt(5.)
    
#    purity = "{:.3f}".format(purity)
#    compl = "{:.3f}".format(compl)
#    merit = "{:.3f}".format(merit)
#    print(f'{class_names[type_orbit]}: {compl}, {purity}, {merit}')

# show the first images in the test data
fig = plt.figure(figsize=(8, 12))
for i in range(50):
    plt.subplot(10, 5, i + 1)
    plt.xticks([])
    plt.yticks([])
    plt.grid(False)
    plt.imshow(test_images[i])
    color = 'blue'
    plt.xlabel("{} ({:2.0f}%)".format(predict_label[i], predict_acc[predict_label[i]]),
               color=color, fontsize=10)
    plt.ylabel("{}".format(test_id[i]), color=color, fontsize=10)
plt.subplots_adjust(hspace=0.3, wspace=0)
# plt.show()
fig.savefig('predicted_data.png', format='png', dpi=300)
plt.close(fig)

print("End")
