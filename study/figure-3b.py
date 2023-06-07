#!/usr/bin/env python3
from __future__ import print_function
import sys
sys.path.append('..')
import os
import numpy as np
import matplotlib.pyplot as plt

import pints.io

import method.heatmap as heatmap  # heap map helper function

import seaborn as sns
cmap = sns.color_palette('flare', as_cmap=True)

#
# Plot all scores for cross criteria for each protocol.
#


# Settings
model_side = ['Single', 'Averaged']
measure_side = ['LSA A', 'LSA D', 'LSA E', 'GSA A', 'GSA D', 'GSA E',]
#                'Shannon']

opt_models = ['ohara-2011', 'model-list']
opt_measures = ['LSA-A', 'LSA-D', 'LSA-E', 'GSA-A', 'GSA-D', 'GSA-E',]
#                'Shannon']

opt_model_labels = ['Single model', 'Average model']


savedir = './fig/'
if not os.path.isdir(savedir):
    os.makedirs(savedir)

inputdir = './practicality-input'

try:
    i_model = int(sys.argv[1])
except IndexError:
    i_model = 3

model_list = ['tnnp', 'fink', 'grandi', 'ohara', 'cipa', 'tomek']

mt = mf = model_list[i_model]
f = '%s/true_%s-fit_%s-row_models-col_measures-vc.txt' \
    % (inputdir, mt, mf)
id_matrix = np.loadtxt(f, dtype=int)

f = '%s/ch3.txt' % inputdir
ch3 = np.loadtxt(f, dtype=int)[i_model]

f = '%s/groenendaal-2015.txt' % inputdir
gro = np.loadtxt(f, dtype=int)[i_model]


# Go through designs
all_std = []
lastniter = 10000
thinning = 2
for ii, model in enumerate(model_side):
    row_std = []
    for jj, measure in enumerate(measure_side):
        int_id = id_matrix[ii, jj]
        run_id = '%03d' % int_id

        # Load samples
        loadas = 'practicality-mcmc/run_%s' % run_id
        s = np.array(pints.io.load_samples('%s-chain.csv' % loadas, n=3))
        s = s[:, -lastniter::thinning, :]
        s = s.reshape(-1, s.shape[-1])
        # std = np.mean(np.std(s, axis=0))
        # Calculat RMSE
        std = np.mean(np.sqrt(np.mean((s - 1)**2, axis=0)))

        row_std.append(std)
    all_std.append(row_std)

benchmark = []
loadas = 'practicality-mcmc/run_%03d' % ch3
s = np.array(pints.io.load_samples('%s-chain.csv' % loadas, n=3))
s = s[:, -lastniter::thinning, :]
s = s.reshape(-1, s.shape[-1])
# std = np.mean(np.std(s, axis=0))
# Calculat RMSE
std = np.mean(np.sqrt(np.mean((s - 1)**2, axis=0)))
benchmark.append(std)

loadas = 'practicality-mcmc/run_%03d' % gro
s = np.array(pints.io.load_samples('%s-chain.csv' % loadas, n=3))
s = s[:, -lastniter::thinning, :]
s = s.reshape(-1, s.shape[-1])
# std = np.mean(np.std(s, axis=0))
# Calculat RMSE
std = np.mean(np.sqrt(np.mean((s - 1)**2, axis=0)))
benchmark.append(std)
del(s)

all_std = np.asarray(all_std) * 1e3
benchmark = np.asarray(benchmark).reshape(1, -1) * 1e3


# Plot heatmap
clim = (np.min(all_std), np.max(all_std))
thres = clim[0] + (clim[1] - clim[0]) * 0.2

fig, axes = plt.subplots(2, 1, gridspec_kw={'height_ratios': [7, 2]},
                               figsize=(6, 4))
ax1, ax2 = axes

im1, _ = heatmap.heatmap(all_std, opt_model_labels, measure_side,
                           ax=ax1, cmap=cmap,
                           #cbarlabel='Averaged ranking score (%)')
                           cbarlabel=None)
_ = heatmap.annotate_heatmap(im1, valfmt='{x:.1f}', threshold=thres)
im1.set_clim(clim)

names = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L']
for i, name in enumerate(names):
    xi, yi = i % 6, int(i / 6)
    x = xi * (1 / 6.) + (1 / 6.) * 0.05
    y = (2 - yi) * (1 / 2.) - (1 / 2.) * 0.05
    color = 'white' if i in [0, 4, 5, 7, 8, 11] else 'black'
    ax1.text(x, y, name, transform=ax1.transAxes,
             ha='left', va='top', weight='bold', color=color)

im2, _ = heatmap.heatmap(benchmark,
                         ['Benchmark'],
                         ['Lei et al.', 'Groenendaal et al.'],
                         ax=ax2, cmap=cmap, cbarlabel=None,
                         rotation=-15)
_ = heatmap.annotate_heatmap(im2, valfmt='{x:.1f}', threshold=thres)
im2.set_clim(clim)

names = ['M', 'N']
for i, name in enumerate(names):
    xi, yi = i % 2, int(i / 2)
    x = xi * (1 / 2.) + (1 / 2.) * 0.05
    y = (1 - yi) * (1 / 1.) - (1 / 1.) * 0.05
    color = 'white'
    ax2.text(x, y, name, transform=ax2.transAxes,
             ha='left', va='top', weight='bold', color=color)

#fig.tight_layout()

#fig.subplots_adjust(right=0.9)
#cbar_ax = fig.add_axes([0.925, 0.15, 0.05, 0.8])
cbar = fig.colorbar(im1, ax=axes.ravel().tolist())
cbar.ax.set_ylabel(r'Averaged posterior RMSE $\times10^3$', rotation=-90, va="bottom")

fig.savefig('%s/fig3b-%s.pdf' % (savedir, mt), format='pdf', bbox_inches='tight')
plt.close('all')
