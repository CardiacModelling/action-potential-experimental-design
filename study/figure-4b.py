#!/usr/bin/env python3
import sys
sys.path.append('..')
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import pints.io
import method.model

#
# Plot the posterior predictive under a simple protocol.
#


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

# Create true model and synthetic data
model_true = method.model.CCModel(info_1.true_model_file, transform=None, dt=dt)
model_true.design([])
times = model_true.times()
parameters_true = np.ones(model_true.n_parameters())
data_no_noise = model_true.simulate(parameters_true)
data = np.copy(data_no_noise)
data += np.random.normal(0, info_1.noise_sigma, size=data.shape)

# Create fitting model
model_pred = method.model.CCModel(info_1.fit_model_file, transform=None, dt=dt)
model_pred.design([])

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


# Plot outputs
fig, ax = plt.subplots(figsize=(4, 2.5))
axins = ax.inset_axes([300., -20., 250, 70], transform=ax.transData)

ax.plot(times, data, c='#7f7f7f', alpha=0.75, label='Data')
axins.plot(times, data, c='#7f7f7f', alpha=0.75)
for i, si in enumerate(s1):
    o = model_pred.simulate(si)
    ax.plot(times, o, alpha=0.25 if i else 0.8, lw=0.5 if i else 1, c='C0',
            label='_' if i else 'Biomarker')
    axins.plot(times, o, alpha=0.25, lw=0.5, c='C0')
for i, si in enumerate(s2):
    o = model_pred.simulate(si)
    ax.plot(times, o, alpha=0.25 if i else 0.8, lw=0.5 if i else 1, c='C1',
            label='_' if i else 'LSA A')
    axins.plot(times, o, alpha=0.25, lw=0.5, c='C1')
ax.plot(times, data_no_noise, alpha=1, lw=0.8, ls='--', c='k',
        label='Groundtruth')
axins.plot(times, data_no_noise, alpha=1, lw=0.8, ls='--', c='k')

ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')

x1, x2 = 40, 110
y1, y2 = 35, 50
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xticks([])
axins.set_yticks([])

rect = patches.Rectangle((x1, y1), x2-x1, y2-y1,
                         lw=0.75, ec='#7f7f7f', fc='none', alpha=0.75)
ax.add_patch(rect)
#ax.indicate_inset_zoom(axins, edgecolor='black')

ax.legend(loc=4, fontsize=8)

plt.tight_layout()

plt.savefig('%s/fig4b.pdf' % (savedir), format='pdf', bbox_inches='tight')
plt.close()
