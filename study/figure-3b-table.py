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

from method.model import parameter_names

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
        std = np.sqrt(np.mean((s - 1)**2, axis=0))

        all_std.append(std)

benchmark = []
loadas = 'practicality-mcmc/run_%03d' % ch3
s = np.array(pints.io.load_samples('%s-chain.csv' % loadas, n=3))
s = s[:, -lastniter::thinning, :]
s = s.reshape(-1, s.shape[-1])
# std = np.mean(np.std(s, axis=0))
# Calculat RMSE
std = np.sqrt(np.mean((s - 1)**2, axis=0))
benchmark.append(std)

loadas = 'practicality-mcmc/run_%03d' % gro
s = np.array(pints.io.load_samples('%s-chain.csv' % loadas, n=3))
s = s[:, -lastniter::thinning, :]
s = s.reshape(-1, s.shape[-1])
# std = np.mean(np.std(s, axis=0))
# Calculat RMSE
std = np.sqrt(np.mean((s - 1)**2, axis=0))
benchmark.append(std)
del(s)

all_std = np.asarray(all_std) * 1e3
benchmark = np.asarray(benchmark) * 1e3


# Make table
#r'Averaged posterior RMSE $\times10^3$'
names = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L']
bench = ['M', 'N']
#_ = heatmap.annotate_heatmap(im2, valfmt='{x:.1f}', threshold=thres)

assert(len(names) == len(all_std))
assert(len(bench) == len(benchmark))

heading = 'Protocol '
for n in parameter_names:
    heading += '& ' + n + ' '
heading += '\\\\'
print(heading)
print('\\midrule')

for n, z in zip(names, all_std):
    l = '\\textbf{' + n +'} '
    for i, s in enumerate(z):
        if i in [0, 6]:
            l += '& %.2f ' % s
        else:
            l += '& %.1f ' % s
    l += '\\\\'
    print(l)

for n, z in zip(bench, benchmark):
    l = '\\textbf{' + n +'} '
    for i, s in enumerate(z):
        if i in [0, 6]:
            l += '& %.2f ' % s
        else:
            l += '& %.1f ' % s
    l += '\\\\'
    print(l)
