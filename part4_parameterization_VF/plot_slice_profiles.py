# %%
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
#%%
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='20')
plt.rc('ytick', labelsize='20')
plt.rc('text', usetex=True)
# %% # --------------------------------------------------------- #
grid = 129
file_dir = os.path.dirname(os.path.realpath(__file__))
case_dir = file_dir + "/W12"
Slice_files = case_dir + "/slice_50"
Slice = np.fromfile(Slice_files, dtype = float).reshape(3, grid, grid)
#%% determine the circle area - create x, y index
xi = np.linspace(0, 100, grid)
yi = np.linspace(0, 100, grid)
by = 30
# Roof parameters (from Green2D.c)
NUM_WAVES = 2
ROOF_BASE = 70.
ROOF_AMPLITUDE = 15.
ROOF_CENTER = ROOF_BASE + ROOF_AMPLITUDE
# Calculate roof shape
X, Y = np.meshgrid(xi, yi)
y_roof = ROOF_CENTER + ROOF_AMPLITUDE * np.cos(2 * np.pi * NUM_WAVES * X / 100.)
# %% plot the cw profile
c_sto_ind = int(grid * 52./100.)
fig1 = plt.figure(figsize=(6, 8))
ax = fig1.add_subplot(1,1,1)
ax.plot(Slice[0,c_sto_ind,:], yi, linewidth = 2)
ax.plot(Slice[0,1,:],  yi, linewidth = 2)
ax.set_xlabel(r'$c_w$', fontsize = 20)
ax.set_ylabel(r'$y$', fontsize = 20)
ax.axhline(y=by, color='k', linestyle='--')
plt.grid(linestyle = ':')
plt.tight_layout()
plt.savefig("profile_cw0.png")
# %% plot the ux profile
fig1 = plt.figure(figsize=(6, 8))
ax = fig1.add_subplot(1,1,1)
ax.plot(Slice[1,c_sto_ind,:], yi, linewidth = 2)
ax.plot(Slice[1,1,:], yi, linewidth = 2)
ax.set_xlabel(r'$u_{x}$', fontsize = 20)
ax.set_ylabel(r'$y$', fontsize = 20)
ax.axhline(y=by, color='k', linestyle='--')
plt.grid(linestyle = ':')
plt.tight_layout()
plt.savefig('profile_Ux.png')
# %% plot 2D slice of cw
fig2 = plt.figure(figsize=(10, 8))
ax = fig2.add_subplot(1,1,1)
# Mask data above roof
cw_data = Slice[0,:,:].T.copy()
cw_data[Y > y_roof] = np.nan
im = ax.contourf(X, Y, cw_data, levels=20, cmap='Blues')
plt.colorbar(im, ax=ax, label=r'$c_w$')
ax.plot(xi, y_roof[0,:], 'k-', linewidth=2, label='Roof')
ax.set_xlabel(r'$x$', fontsize = 20)
ax.set_ylabel(r'$y$', fontsize = 20)
ax.axhline(y=by, color='k', linestyle='--')
plt.tight_layout()
plt.savefig('2D_slice_cw.png')
# %% plot 2D slice of ux
fig3 = plt.figure(figsize=(10, 8))
ax = fig3.add_subplot(1,1,1)
# Mask data above roof
ux_data = Slice[1,:,:].T.copy()
ux_data[Y > y_roof] = np.nan
im = ax.contourf(X, Y, ux_data, levels=20, cmap='RdBu_r')
plt.colorbar(im, ax=ax, label=r'$u_x$')
ax.plot(xi, y_roof[0,:], 'k-', linewidth=2, label='Roof')
ax.set_xlabel(r'$x$', fontsize = 20)
ax.set_ylabel(r'$y$', fontsize = 20)
ax.axhline(y=by, color='k', linestyle='--')
plt.tight_layout()
plt.savefig('2D_slice_ux.png')
# %%
