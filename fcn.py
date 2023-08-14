import tensorflow as tf
import numpy as np
import h5py
import os

# Parameters
data_path = r'C:\Users\ofirl\PycharmProjects\pythonProject\pythonProject\NCBI\e_coli'
data_path1 = r'C:\Users\ofirl\PycharmProjects\pythonProject\pythonProject\NCBI\staphylococcus_aureus'
input_size = 10  # Specify the input size based on your data
output_size = 2  # Specify the output size based on the number of classes
batch_size = 32

READ_SIZE = 1000000
NUMBER_OF_READS = 1

# Load data from HDF files
x_data = []
y_data = []

x_test = []
y_test = []

flag = 0
for fn in os.listdir(data_path):
    if os.path.isfile(os.path.join(data_path, fn)) and (os.path.splitext(fn)[-1].lower() == ".hdf5"):
        with h5py.File(os.path.join(data_path, fn), 'r') as f:
            read = []
            for i in range(NUMBER_OF_READS):
                s = "read_" + str(i)
                # read.append(f[s][()])
                read = f[s][()]
                if flag < 1:
                    x_test.append(read)
                    y_test.append(0)
                    flag = flag + 1
                else:
                    x_data.append(read)
                    y_data.append(0)
            # read = np.asarray(read)
            # read = read.flatten()
            # read = tf.convert_to_tensor(read, dtype=tf.float32)
            # if flag < 1:
            #     x_test.append(read)
            #     y_test.append(0)
            #     flag = flag + 1
            # else:
            #     x_data.append(read)
            #     y_data.append(0)

flag = 0
for fn in os.listdir(data_path1):
    if os.path.isfile(os.path.join(data_path1, fn)) and (os.path.splitext(fn)[-1].lower() == ".hdf5"):
        with h5py.File(os.path.join(data_path1, fn), 'r') as f:
            read = []
            for i in range(NUMBER_OF_READS):
                s = "read_" + str(i)
                # read.append(f[s][()])
                read = f[s][()]
                if flag < 1:
                    x_test.append(read)
                    y_test.append(1)
                    flag = flag + 1
                else:
                    x_data.append(read)
                    y_data.append(1)
            # read = np.asarray(read)
            # read = read.flatten()
            # read = tf.convert_to_tensor(read, dtype=tf.float32)
            # if flag < 1:
            #     x_test.append(read)
            #     y_test.append(1)
            #     flag = flag + 1
            # else:
            #     x_data.append(read)
            #     y_data.append(1)

# randomise samples
n_samples = len(y_data)
rand_gen = np.random.RandomState(0)
indices = np.arange(n_samples)
rand_gen.shuffle(indices)
x_data = [x_data[i] for i in indices]
y_data = [y_data[i] for i in indices]

# Convert data to NumPy arrays
x_data = np.asarray(x_data)
y_data = np.asarray(y_data)
x_test = np.asarray(x_test)
y_test = np.asarray(y_test)

# x_data = np.expand_dims(x_data, axis=1)
# x_test = np.expand_dims(x_test, axis=1)
y_data = np.expand_dims(y_data, axis=1)
y_test = np.expand_dims(y_test, axis=1)

y_train_binary = tf.keras.utils.to_categorical(y_data, output_size)
y_test_binary = tf.keras.utils.to_categorical(y_test, output_size)

# Create a fully connected neural network model
model = tf.keras.Sequential([
    tf.keras.layers.Dense(units=512, activation='relu', input_shape=x_data[0].shape),
    tf.keras.layers.Dense(units=128, activation='relu'),
    tf.keras.layers.Dense(units=output_size, activation='softmax')
])

# Compile the model
model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])

# Train the model
model.fit(x_data, y_train_binary, epochs=10, batch_size=2)

predictions = model.predict(x_test)

[print(*line) for line in predictions]

# Save the trained model
model.save('trained_model.h5')
