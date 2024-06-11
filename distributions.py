import numpy as np
import scipy
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from tqdm import tqdm

ALLrs = np.zeros(0)
ALLangles = np.zeros(0)

def do(As, Bs, Cs, Ds):
    Qs = 0.5*(Cs+Ds)
    Ws = 0.5*(As+Bs)

    RQWs = np.sqrt(np.sum(np.square(Qs-Ws), axis=1))

    RABs = np.sqrt(np.sum(np.square(As-Bs), axis=1))
    RCDs = np.sqrt(np.sum(np.square(Cs-Ds), axis=1))
    angles = np.arccos(np.sum((As-Bs)*(Cs-Ds), axis=1)/(RABs*RCDs))


    RL = np.where(RABs > RCDs, RABs, RCDs)
    mask = (RQWs < 0.01) & (RL > 0.1) & (RL < 0.2)

    rs = np.where(RABs < RCDs, RABs/RCDs, RCDs/RABs)
    #rs = np.where(RABs < RCDs, RCDs, RABs)

    global ALLrs, ALLangles
    ALLrs = np.concatenate([ALLrs, rs[mask]])
    ALLangles = np.concatenate([ALLangles, angles[mask]])

for REP in tqdm(range(1000)):
    N = 100000

    dim = 2

    As = np.random.normal(0, 0.4, size=(N, dim))
    Bs = np.random.normal(0, 0.4, size=(N, dim))
    Cs = np.random.normal(0, 0.4, size=(N, dim))
    Ds = np.random.normal(0, 0.4, size=(N, dim))

    do(As, Bs, Cs, Ds)
    do(As, Cs, Bs, Ds)
    do(As, Ds, Bs, Cs)

cmap = 'viridis'
norm = LogNorm

fig,ax = plt.subplots(figsize=(20,20))
plt.hist(ALLrs, bins=20)
plt.show()
fig,ax = plt.subplots(figsize=(20,20))
plt.hist(ALLangles, bins=20)
plt.show()
fig,ax = plt.subplots(figsize=(20,20))
plt.hist2d(ALLangles, ALLrs, bins=20,
           cmap=cmap, norm=norm())
plt.show()
fig,ax = plt.subplots(figsize=(20,20), subplot_kw={'projection': 'polar'})
ax.hist2d(ALLangles, ALLrs, bins=20,
          cmap=cmap, norm=norm())
plt.show()
