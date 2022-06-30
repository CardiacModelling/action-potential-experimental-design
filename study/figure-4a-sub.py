#!/usr/bin/env python3
import sys
sys.path.append('..')
import os
import numpy as np
import matplotlib.pyplot as plt

import pints.io
import pints.plot
import method.model

savedir = './fig/'
if not os.path.isdir(savedir):
    os.makedirs(savedir)

input_file = './practicality-input/id_160.py'

inputdir = './practicality-input'

# Get all input variables
import importlib
sys.path.append(os.path.dirname(input_file))
base = os.path.basename(input_file)
info = importlib.import_module(os.path.splitext(base)[0])

print('Run ID: ', info.run_id)

# Settings
dt = 0.1  # ms
seed_id = 101  # random seed
np.random.seed(seed_id)
fit_seed = np.random.randint(0, 2**30)
np.random.seed(fit_seed)
n_steps = 0

# Create true model and synthetic data
model_true = method.model.CCBiomarkerModel(info.true_model_file, transform=None, dt=dt)
times = model_true.times()
parameters_true = np.ones(model_true.n_parameters())
raw_data_no_noise = model_true._simulate(parameters_true)
model_true._biomarkers.set_data(model_true._pacing, times, raw_data_no_noise)
data_no_noise = model_true.extract()

data_list = []
np.random.seed(0)
for i in range(200):
    raw_data = np.copy(raw_data_no_noise)
    raw_data += np.random.normal(0, info.noise_sigma, size=raw_data.shape)
    model_true._biomarkers.set_data(model_true._pacing, times, raw_data)
    data = model_true.extract()
    data_list.append(data)
data_list = np.array(data_list)

# Create fitting model
# model_pred = method.model.CCBiomarkerModel(info.fit_model_file, transform=None, dt=dt)


fig, axes = plt.subplots(1, 4, figsize=(8, 1.75), sharey=True)
axes = axes.reshape(1, -1)
nx, ny = axes.shape
for j, k in enumerate(['mur', 'apa', 'vmax', 'apd10']):
    a, b = int(j / ny), j % ny
    i = model_true._biomarkers.list_of_biomarkers.index(k)
    xmin, xmax = np.min(data_list[:, i]), np.max(data_list[:, i])
    xmin = min(xmin, data_no_noise[i])
    xmax = max(xmax, data_no_noise[i])
    axes[a, b].hist(data_list[:, i], bins=np.linspace(xmin, xmax, 30))
    axes[a, b].axvline(data_no_noise[i], c='#7f7f7f', ls='--')
    axes[a, b].set_xlabel(model_true._biomarkers.biomarker_short_name(k)
                          + ' ('
                          + model_true._biomarkers.biomarker_unit(k)
                          + ')')
    if b == 0:
        axes[a, b].set_ylabel('Frequency')

#axes[-1, -1].axis('off')

plt.tight_layout()

plt.savefig('%s/fig4a-sub.pdf' % (savedir), format='pdf', bbox_inches='tight')
plt.close()
