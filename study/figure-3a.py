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


to_plot = [3, 1, -2, -1]  # best two and worst two
names = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']

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
all_samples = []
all_names = []
lastniter = 10000
thinning = 2
for ii, model in enumerate(model_side):
    for jj, measure in enumerate(measure_side):
        int_id = id_matrix[ii, jj]
        run_id = '%03d' % int_id

        # Load samples
        loadas = 'practicality-mcmc/run_%s' % run_id
        s = np.array(pints.io.load_samples('%s-chain.csv' % loadas, n=3))
        s = s[:, -lastniter::thinning, :]
        s = s.reshape(-1, s.shape[-1])

        all_samples.append(s)
        all_names.append(model + ', ' + measure)
        del(s)

loadas = 'practicality-mcmc/run_%03d' % ch3
s = np.array(pints.io.load_samples('%s-chain.csv' % loadas, n=3))
s = s[:, -lastniter::thinning, :]
s = s.reshape(-1, s.shape[-1])

all_samples.append(s)
all_names.append('Lei et al.')
del(s)

loadas = 'practicality-mcmc/run_%03d' % gro
s = np.array(pints.io.load_samples('%s-chain.csv' % loadas, n=3))
s = s[:, -lastniter::thinning, :]
s = s.reshape(-1, s.shape[-1])

all_samples.append(s)
all_names.append('Groenendaal et al.')
del(s)


all_samples = [all_samples[i] for i in to_plot]
#all_names = [all_names[i] for i in to_plot]
all_names = ['Protocol ' + names[i] for i in to_plot]


# Plot histograms
bins = 40
alpha = 0.75
n_percentiles = None

fig, axes = plt.subplots(2, 4, figsize=(10, 4))
plt.subplots_adjust(hspace=.45, wspace=.25)
for i in range(axes.size):
    ai, aj = int(i // 4), i % 4

    axes[ai, aj].axvline(1, color='k', linestyle='--', label='Ground truth')

    axes[ai, aj].set_xlabel(parameters_nice[i], fontsize=14)
    if aj == 0:
        axes[ai, aj].set_ylabel('Normalised\nposterior', fontsize=12)

    axes[ai, aj].ticklabel_format(axis='both', style='sci', scilimits=(-2, 3))

    for j, samples_j in enumerate(all_samples):
        if n_percentiles is None:
            xmin = np.min(samples_j[:, i])
            xmax = np.max(samples_j[:, i])
        else:
            xmin = np.percentile(samples_j[:, i],
                                 50 - n_percentiles / 2.)
            xmax = np.percentile(samples_j[:, i],
                                 50 + n_percentiles / 2.)
        xbins = np.linspace(xmin, xmax, bins)

        #axes[ai, aj].hist(samples_j[:, i], bins=xbins,
        #        alpha=alpha, histtype='step', linewidth=1.5,
        #        density=True, label=all_names[j], color='C' + str(j),
        #        zorder=-j)
        H, _ = np.histogram(samples_j[:, i], bins=xbins)
        axes[ai, aj].plot(xbins[:-1], H/np.max(H), ds='steps', label=all_names[j], color='C' + str(j), alpha=alpha, linewidth=1.5, zorder=-j)

    axes[ai, aj].set_ylim([0, 1.02])

axes[0, 0].legend(loc='lower left', bbox_to_anchor=(-0.05, 1.1), ncol=5,
        bbox_transform=axes[0, 0].transAxes)
plt.savefig('%s/fig3a-%s.pdf' % (savedir, mt), format='pdf', bbox_inches='tight')
plt.close()
