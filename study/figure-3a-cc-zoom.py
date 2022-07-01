#!/usr/bin/env python3
from __future__ import print_function
import sys
sys.path.append('..')
import os
import numpy as np
import matplotlib.pyplot as plt

import pints.io

from method.model import parameter_names as parameters_nice

#
# Plot all scores for cross criteria for each protocol.
#


to_plot = [10, 9, -2, -1]  # best two and worst two
names = ['U', 'V', 'M', 'N', 'O', 'P', 'Q', 'R', 'U', 'V']
skip_chain = {  # ID, chain ID
    '123': [1],
    #'166': [0],
}

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
bm = '%03d' % np.loadtxt(f, dtype=int)[3]

f = '%s/1hz.txt' % inputdir
hz1 = '%03d' % np.loadtxt(f, dtype=int)[3]


# Go through designs
all_samples = []
all_names = []
lastniter = 10000
thinning = 2

loadas = 'practicality-mcmc-bm/run_%s' % bm
s = np.array(pints.io.load_samples('%s-chain.csv' % loadas, n=3))
s = s[:, -lastniter::thinning, :]
s = s.reshape(-1, s.shape[-1])

all_samples.append(s)
all_names.append('Biomarkers')
del(s)

loadas = 'practicality-mcmc-cc/run_%s' % hz1
s = np.array(pints.io.load_samples('%s-chain.csv' % loadas, n=3))
print( hz1, hz1 in skip_chain.keys())
if hz1 in skip_chain.keys():
    s = np.delete(s, skip_chain[hz1], axis=0)
s = s[:, -lastniter::thinning, :]
s = s.reshape(-1, s.shape[-1])

all_samples.append(s)
all_names.append('1 Hz')
del(s)

for ii, model in enumerate(model_side):
    for jj, measure in enumerate(measure_side):
        int_id = id_matrix[ii, jj]
        run_id = '%03d' % int_id

        # Load samples
        loadas = 'practicality-mcmc-cc/run_%s' % run_id
        s = np.array(pints.io.load_samples('%s-chain.csv' % loadas, n=3))
        if run_id in skip_chain.keys():
            s = np.delete(s, skip_chain[run_id], axis=0)
        s = s[:, -lastniter::thinning, :]
        s = s.reshape(-1, s.shape[-1])

        all_samples.append(s)
        all_names.append(model + ', ' + measure)
        del(s)


#all_samples = [all_samples[i] for i in to_plot]
#all_names = [all_names[i] for i in to_plot]
all_names = ['Protocol ' + name for name in names]


# Plot histograms
bins = 60
alpha = 0.75
n_percentiles = None

fig, axes = plt.subplots(2, 4, figsize=(10, 4))
plt.subplots_adjust(hspace=.45, wspace=.25)
for i in range(axes.size):
    ai, aj = int(i // 4), i % 4

    axes[ai, aj].axvline(1, color='k', linestyle='--', label='Ground truth',
                         zorder=10)

    axes[ai, aj].set_xlabel(parameters_nice[i], fontsize=14)
    if aj == 0:
        axes[ai, aj].set_ylabel('Marginal\nposterior', fontsize=14)

    axes[ai, aj].ticklabel_format(axis='both', style='sci', scilimits=(-2, 3))

    for j, samples_j in enumerate(all_samples):
        xmin, xmax = 0.98, 1.02

        xbins = np.linspace(xmin, xmax, bins)

        #axes[ai, aj].hist(samples_j[:, i], bins=xbins, alpha=alpha, histtype='step', linewidth=1.5,
        #        density=True, label=all_names[j], color='C' + str(j))
        H, _ = np.histogram(samples_j[:, i], bins=xbins)
        if j in [0, 1]:
            axes[ai, aj].bar(xbins[:-1], H/np.max(H), width=(xbins[1]-xbins[0]), label=all_names[j], color='C' + str(j), alpha=0.5)
        else:
            axes[ai, aj].plot(xbins[:-1], H/np.max(H), ds='steps', label=all_names[j], color='C' + str(j))


axes[0, 0].legend(loc='lower left', bbox_to_anchor=(-0.2, 1.1), ncol=5,
        bbox_transform=axes[0, 0].transAxes)
plt.savefig('%s/fig3a-cc-zoom.pdf' % (savedir), format='pdf', bbox_inches='tight')
plt.close()
