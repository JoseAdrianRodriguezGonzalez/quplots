import json
import numpy as np


def save_wavefunction(x, y, z, psi, filepath):
    data = {
        "x": np.asarray(x).tolist(),
        "y": np.asarray(y).tolist(),
        "z": np.asarray(z).tolist(),
        "psi": np.asarray(psi).tolist(),
    }
    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(data, f)
