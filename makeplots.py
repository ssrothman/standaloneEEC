import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

# Load the data
with open("resolved4_shapes.dat", 'r') as f:
    l0 = f.readline()
    l1 = f.readline()

# Parse the data
l1 = l1.split(', ')
l1 = list(map(float, l1))

DATA = np.asarray(l1)

# Parse the shape
l0 = l0[1:-2]
l0 = l0.split(', ')
SHAPE = tuple(map(int, l0))

DATA = np.reshape(DATA, SHAPE)

print(DATA.shape)

titles = ['No shape', 'Dipole', 'Tee', 'Triangle']

nbin_r = SHAPE[2]
nbin_phi = SHAPE[3]
edges_r = np.linspace(0, 1, nbin_r - 1)
edges_r = np.concatenate(([-np.inf], edges_r, [np.inf]))
edges_phi = np.linspace(0, 0.5*np.pi, nbin_phi - 1)
edges_phi = np.concatenate(([-np.inf], edges_phi, [np.inf]))
centers_r = 0.5*(edges_r[1:] + edges_r[:-1])
centers_phi = 0.5*(edges_phi[1:] + edges_phi[:-1])
print(SHAPE)
print(edges_r)
print(edges_phi)

# Plot the data
for shapenum in [1, 2]:
    plt.figure(figsize=(20, 20))
    plt.title(titles[shapenum], fontsize=40)
    plt.imshow(DATA[shapenum, 2][1:-1, 1:-1],
               norm=Normalize(), origin='lower')
    plt.colorbar()
    plt.xlabel("phi bin", fontsize=30)
    plt.ylabel("r bin", fontsize=30)
    plt.show()
