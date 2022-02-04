#!/usr/bin/env python3
from __future__ import print_function
import sys
sys.path.append('..')
import os
import numpy as np
import matplotlib.pyplot as plt

import pints.io

import method.heatmap as heatmap  # heap map helper function

#
# Plot all scores for cross criteria for each protocol.
#


# Settings
model_side = ['Single', 'Averaged']
measure_side = ['LSA A', 'LSA D', 'LSA E']
#, 'GSA A', 'GSA D', 'GSA E',]
#                'Shannon']

opt_models = ['ohara-2011', 'model-list']
opt_measures = ['LSA-A', 'LSA-D', 'LSA-E']
#, 'GSA-A', 'GSA-D', 'GSA-E',]
#                'Shannon']

opt_model_labels = ['Single model', 'Average model']


savedir = './fig/'
if not os.path.isdir(savedir):
    os.makedirs(savedir)

inputdir = './practicality-input'

mt = mf = 'ohara'
f = '%s/true_%s-fit_%s-row_models-col_measures-cc.txt' \
    % (inputdir, mt, mf)
id_matrix = np.loadtxt(f, dtype=int)

f = '%s/biomarkers.txt' % inputdir
bm = np.loadtxt(f, dtype=int)[3]


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
        loadas = 'practicality-mcmc-cc/run_%s' % run_id
        s = np.array(pints.io.load_samples('%s-chain.csv' % loadas, n=3))
        s = s[:, -lastniter::thinning, :]
        s = s.reshape(-1, s.shape[-1])
        std = np.mean(np.std(s, axis=0))

        row_std.append(std)
    all_std.append(row_std)

benchmark = []
loadas = 'practicality-mcmc-bm/run_%03d' % bm
s = np.array(pints.io.load_samples('%s-chain.csv' % loadas, n=3))
s = s[:, -lastniter::thinning, :]
s = s.reshape(-1, s.shape[-1])
std = np.mean(np.std(s, axis=0))
benchmark.append(std)

del(s)

all_std = np.asarray(all_std) * 1e3
benchmark = np.asarray(benchmark).reshape(1, -1) * 1e3


# Plot heatmap
clim = (np.min(all_std), np.max(all_std))
thres = clim[0] + (clim[1] - clim[0]) * 0.6

fig, axes = plt.subplots(2, 1, gridspec_kw={'height_ratios': [7, 2]},
                               figsize=(6, 4))
ax1, ax2 = axes

im1, _ = heatmap.heatmap(all_std, opt_model_labels, measure_side,
                           ax=ax1, cmap='YlGn',
                           #cbarlabel='Averaged ranking score (%)')
                           cbarlabel=None)
_ = heatmap.annotate_heatmap(im1, valfmt='{x:.1f}', threshold=thres)
im1.set_clim(clim)

im2, _ = heatmap.heatmap(benchmark,
                         ['Benchmark'],
                         ['Biomarkers (VC only)'],
                         ax=ax2, cmap='YlGn', cbarlabel=None,
                         rotation=-15)
_ = heatmap.annotate_heatmap(im2, valfmt='{x:.1f}', threshold=thres)
im2.set_clim(clim)

#fig.tight_layout()

#fig.subplots_adjust(right=0.9)
#cbar_ax = fig.add_axes([0.925, 0.15, 0.05, 0.8])
cbar = fig.colorbar(im1, ax=axes.ravel().tolist())
cbar.ax.set_ylabel(r'Averaged posterior standard deviation $\times10^3$', rotation=-90, va="bottom")

fig.savefig('%s/fig3b-cc.pdf' % (savedir), format='pdf', bbox_inches='tight')
plt.close('all')
