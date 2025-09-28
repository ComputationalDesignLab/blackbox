import numpy as np

data = np.loadtxt("force_field.txt", skiprows=1)  # skip the 1309 1008 line
coords = data[:, :3]
forces = data[:, 3:]

# unique to tolerance (say 1e-6)
# rounded = np.round(coords, 10)
uniq_coords, idx = np.unique(coords, axis=0, return_index=True)
print("rows in file:", coords.shape[0])
print("unique coords (rounded):", uniq_coords.shape[0])
