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
Slice_files = case_dir + "/slice_100"
Slice = np.fromfile(Slice_files, dtype = float).reshape(3, grid, grid)
#%% determine the circle area - create x, y index
xi = np.linspace(0, 100, grid)
yi = np.linspace(0, 100, grid)
by = 30
# %% plot the cw profile
c_sto_ind = int(grid * 0.35)
fig1 = plt.figure(figsize=(6, 8))
ax = fig1.add_subplot(1,1,1)
ax.plot(Slice[0,c_sto_ind,:], yi, linewidth = 2)
ax.plot(Slice[0,1,:],  yi, linewidth = 2)
ax.set_xlabel(r'$c_w$', fontsize = 20)
ax.set_ylabel(r'$y$', fontsize = 20)
ax.axhline(y=by, color='k', linestyle='--')
plt.grid(linestyle = ':')
plt.tight_layout()
# plt.savefig(dir_out + "profile_cw0.png")
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
# plt.savefig(dir_out + 'profile_Ux.png')
# %%
