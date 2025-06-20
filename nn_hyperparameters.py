nn_topology_configs = {
        "GTR": {
           "n_layers": 4,
           "layer_params": [
               {"n_units": 480, "activation": "relu", "dropout":0.17669514670422368},
               {"n_units": 352, "activation": "tanh", "dropout":0.18464804508715815},
               {"n_units": 96, "activation": "sigmoid", "dropout":0.19029944798800483},
               {"n_units": 512, "activation": "elu", "dropout":0.3021441985631763},
           ],
            "optimizer": "adamax",
            "batch_size": 128,
            "learning_rate": 0.00026900530606856465,
       },
}
