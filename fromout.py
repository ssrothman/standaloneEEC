import numpy as np
from matplotlib.colors import Normalize, LogNorm
import matplotlib.pyplot as plt

# Load the data
with open("out", 'r') as f:
    lines = f.readlines()

# Parse the data
thetas = np.asarray([float(line.strip().split(', ')[0]) for line in lines])
rs = np.asarray([float(line.strip().split(', ')[1]) for line in lines])

thetabins = np.asarray([int(line.strip().split(', ')[2]) for line in lines])
rbins = np.asarray([int(line.strip().split(', ')[3]) for line in lines])

mask = rs < 1
thetas = thetas[mask]
rs = rs[mask]
thetabins = thetabins[mask]
rbins = rbins[mask]

cmap = 'viridis'
norm = Normalize

plt.subplots(figsize=(20,20), subplot_kw={'projection': 'polar'})
allthetas = np.concatenate([thetas, np.pi-thetas, thetas+np.pi, 2*np.pi-thetas])
allrs = np.concatenate([rs, rs, rs, rs])
plt.hist2d(allthetas, allrs, bins=30,
           cmap=cmap, norm=norm())
plt.colorbar()
plt.show()

plt.subplots(figsize=(20,20))
plt.hist2d(thetabins, rbins, bins=[np.max(thetabins),np.max(rbins)],
           cmap=cmap, norm=norm())
plt.colorbar()
plt.show()
