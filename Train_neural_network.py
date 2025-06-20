import sys
import numpy as np
import tensorflow as tf
from sklearn.model_selection import train_test_split
from nn_hyperparameters import nn_topology_configs

# Command line arguments
path_to_frequencies = sys.argv[1]
path_to_labels = sys.argv[2]

# Load training data
frequency_array = np.load(path_to_frequencies)
topology_array = np.load(path_to_labels)

# Preparation (split) of the training data
frequency_train, frequency_test, topology_train, topology_test = train_test_split(
    frequency_array, topology_array, test_size=0.04, random_state=42
)

frequency_train = np.asarray(frequency_train)
topology_train = np.asarray(topology_train)

frequency_test = np.asarray(frequency_test)
topology_test = np.asarray(topology_test)

topology_train = tf.keras.utils.to_categorical(topology_train, 3)
topology_test = tf.keras.utils.to_categorical(topology_test, 3)

# Training data normalization
normlayer = tf.keras.layers.experimental.preprocessing.Normalization()
normlayer.adapt(frequency_train)


# Function to build neural network based on the file with hyperparameters
def build_nn_model(config, normlayer):
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.InputLayer(input_shape=256, ))
    model.add(normlayer)

    for layer_param in config["layer_params"]:
        model.add(
            tf.keras.layers.Dense(
                units=layer_param["n_units"],
                activation=layer_param["activation"]
            )
        )
        model.add(tf.keras.layers.Dropout(layer_param["dropout"]))

    model.add(tf.keras.layers.Dense(units=3, activation="softmax"))

    # Extract learning rate if present
    learning_rate = config.get("learning_rate", None)

    if learning_rate:
        optimizer = getattr(tf.keras.optimizers, config["optimizer"].capitalize())(learning_rate=learning_rate)
    else:
        optimizer = getattr(tf.keras.optimizers, config["optimizer"].capitalize())()

    model.compile(optimizer=optimizer, loss="categorical_crossentropy", metrics=["accuracy"])
    return model


def train_model(
    model_func, X_train, Y_train, X_test, Y_test, normlayer):

    config = nn_topology_configs["GTR"]

    model = model_func(config,normlayer)

    batch_size = nn_topology_configs["GTR"]["batch_size"]

    result = model.fit(
            x=X_train,
            y=Y_train,
            batch_size=batch_size,
            epochs=5,
            validation_data=(X_test, Y_test),
            verbose=2,
        )


    return model



fitted_model = train_model(build_nn_model, frequency_train, topology_train, frequency_test, topology_test, normlayer)

fitted_model.save("MySuperCoolModel")