import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

# Load the data
with open("dipole.dat", 'r') as f:
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

nbin_r = SHAPE[1]
nbin_phi = SHAPE[2]
edges_r = np.linspace(0, 1, nbin_r - 1)
edges_phi = np.linspace(0, 0.5*np.pi, nbin_phi - 1)
centers_r = 0.5*(edges_r[1:] + edges_r[:-1])
centers_phi = 0.5*(edges_phi[1:] + edges_phi[:-1])
print(SHAPE)
print(edges_r)
print(edges_phi)

# Plot the data
fig, ax = plt.subplots(figsize=(20,20), 
                       subplot_kw={'projection': 'polar'})

plt.title("DIPOLE", fontsize=40)

#normclass = LogNorm
normclass = Normalize
cmap = 'viridis'

DATA = DATA[1]

#area = area * np.exp(-np.square(centers_r[:,None])/0.4)
area = np.ones_like(DATA[1:-1, 1:-1])
area = np.pi * (edges_r[1:,None]**2 - edges_r[:-1,None]**2) * (edges_phi[None,1:] - edges_phi[None,:-1]) / (2*np.pi)

DATA = DATA[1:-1, 1:-1]/area
DATA = DATA/np.mean(DATA)

pc = ax.pcolormesh(edges_phi, edges_r, 
                   DATA,
                   cmap=cmap, norm=Normalize(vmin=0.8, vmax=1.2))
pc2 = ax.pcolormesh(np.pi-edges_phi, edges_r, 
                    DATA,
                   cmap=cmap, norm=Normalize(vmin=0.8, vmax=1.2))
pc3 = ax.pcolormesh(np.pi+edges_phi, edges_r, 
                   DATA,
                   cmap=cmap, norm=Normalize(vmin=0.8, vmax=1.2))
pc4 = ax.pcolormesh(2*np.pi-edges_phi, edges_r, 
                    DATA,
                   cmap=cmap, norm=Normalize(vmin=0.8, vmax=1.2))
fig.colorbar(pc)

#plt.imshow(DATA[shapenum, 2][1:-1, 1:-1],
#           norm=Normalize(), origin='lower')
#plt.colorbar()
#plt.xlabel("phi bin", fontsize=30)
#plt.ylabel("r bin", fontsize=30)
plt.show()

#
