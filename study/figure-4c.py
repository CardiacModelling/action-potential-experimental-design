#!/usr/bin/env python3
from __future__ import print_function
import sys
sys.path.append('..')
import os
import numpy as np
import matplotlib.pyplot as plt

import pints.io
import method.model

#
# Plot the posterior predictive under a simple protocol.
#

block = np.zeros(8)
block[2] = 0.9  # ikr
block[3] = 0.55  # iks
block = 1 - block

skip_chain = {'123': [1]}  # ID, chain ID

#try:
#    input_file = sys.argv[1]
#except:
#    print('Usage: python %s [str:input_file]' % os.path.basename(__file__))
#    sys.exit()
input_file_1 = './practicality-input/id_160.py'
input_file_2 = './practicality-input/id_121.py'

savedir = './fig/'
if not os.path.isdir(savedir):
    os.makedirs(savedir)

inputdir = './practicality-input'

# Get all input variables
import importlib
sys.path.append(os.path.dirname(input_file_1))
base = os.path.basename(input_file_1)
info_1 = importlib.import_module(os.path.splitext(base)[0])

sys.path.append(os.path.dirname(input_file_2))
base = os.path.basename(input_file_2)
info_2 = importlib.import_module(os.path.splitext(base)[0])

print('Run ID 1: ', info_1.run_id)
print('Run ID 2: ', info_2.run_id)

# Settings
dt = 0.1  # ms
seed_id = 101  # random seed
np.random.seed(seed_id)
fit_seed = np.random.randint(0, 2**30)
np.random.seed(fit_seed)
n_steps = 0

pacing = [50, 1050, 2050, 3050]

# Create true model and synthetic data
model_true = method.model.CCModel(info_1.true_model_file, transform=None, dt=dt)
model_true.design(pacing)
times = model_true.times()
parameters_true = np.ones(model_true.n_parameters()) * block
data_no_noise = model_true.simulate(parameters_true)
data = np.copy(data_no_noise)
data += np.random.normal(0, info_1.noise_sigma, size=data.shape)

# Create fitting model
model_pred = method.model.CCModel(info_1.fit_model_file, transform=None, dt=dt)
model_pred.design(pacing)

# NOTE: Assume both pred models are the same.

# Go through designs
lastniter = 2000
thinning = 100

# Load samples
try:
    loadas = 'practicality-mcmc-cc/run_%s' % info_1.run_id
    s1 = np.array(pints.io.load_samples('%s-chain.csv' % loadas, n=3))
except:
    loadas = 'practicality-mcmc-bm/run_%s' % info_1.run_id
    s1 = np.array(pints.io.load_samples('%s-chain.csv' % loadas, n=3))
if info_1.run_id in skip_chain.keys():
    s1 = np.delete(s1, skip_chain[info_1.run_id], axis=0)
s1 = s1[:, -lastniter::thinning, :]
s1 = s1.reshape(-1, s1.shape[-1])
print(s1.shape)

try:
    loadas = 'practicality-mcmc-cc/run_%s' % info_2.run_id
    s2 = np.array(pints.io.load_samples('%s-chain.csv' % loadas, n=3))
except:
    loadas = 'practicality-mcmc-bm/run_%s' % info_2.run_id
    s2 = np.array(pints.io.load_samples('%s-chain.csv' % loadas, n=3))
if info_2.run_id in skip_chain.keys():
    s2 = np.delete(s2, skip_chain[info_2.run_id], axis=0)
s2 = s2[:, -lastniter::thinning, :]
s2 = s2.reshape(-1, s2.shape[-1])
print(s2.shape)

r = times < 6000

# Plot histograms
fig,ax = plt.subplots(figsize=(4, 2.5))
ax.plot(times[r], data[r], c='#7f7f7f', alpha=0.75, label='Data')
for i, si in enumerate(s1):
    o = model_pred.simulate(si * block)
    ax.plot(times[r], o[r],
            alpha=0.25 if i else 0.8, lw=0.5 if i else 1, c='C0',
            label='_' if i else 'Biomarker')
for i, si in enumerate(s2):
    o = model_pred.simulate(si * block)
    ax.plot(times[r], o[r],
            alpha=0.25 if i else 0.8, lw=0.5 if i else 1, c='C1',
            label='_' if i else 'LSA A')
ax.plot(times[r], data_no_noise[r], alpha=1, lw=0.8, ls='--', c='k',
        label='Groundtruth')

ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')

plt.tight_layout()

plt.savefig('%s/fig4c.pdf' % (savedir), format='pdf', bbox_inches='tight')
plt.close()
