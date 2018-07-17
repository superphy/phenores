import pickle
import numpy as np
from keras.layers import Dense
from keras.models import Sequential
from keras.utils import to_categorical
from keras.wrappers.scikit_learn import KerasClassifier

with open(snakemake.input[0], 'rb') as f:
    train_data = pickle.load(f)
with open(snakemake.input[1], 'rb') as f:
    train_labels = pickle.load(f)
with open(snakemake.input[2], 'rb') as f:
    test_labels = pickle.load(f)

all_files = np.concatenate((train_labels, test_labels), axis=0)
num_classes = np.unique(all_files).shape[0]

train_labels = to_categorical(train_labels)

input_shape = train_data.shape[1]

model = Sequential()
model.add(Dense(input_shape, input_dim=input_shape, kernel_initializer='normal',
                activation='relu'))
model.add(Dense(num_classes, kernel_initializer='normal'))
model.compile(loss='mean_squared_error', optimizer='adam', metrics=['accuracy'])

model.fit(train_data, train_labels, epochs=100, batch_size=10, verbose=1)

model.save(snakemake.output[0])
