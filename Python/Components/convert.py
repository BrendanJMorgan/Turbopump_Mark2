import numpy as np
import os
path = os.path.dirname(__file__)
# Load the npz file
data = np.load(os.path.join(path, 'impulse_turbine_efficiency.npz'))

print(data.files)   # See what key(s) exist

# Assume the array is stored under the first key
arr = data[data.files[0]]

# Save it as .npy
np.save(os.path.join(path, 'impulse_turbine_efficiency.npy'), arr)
