#!/usr/bin/env python3
"""
# Plotting voltage-clamp OED results.
"""
import sys
sys.path.append('..')
import os
import numpy as np
import matplotlib.pyplot as plt

import method.model


# Settings
n_steps = 20  # number of steps of the protocol
dt = 5  # ms

model_side = ['Single model', 'Averaged model']
measure_side = ['LSA A', 'LSA D', 'LSA E', 'GSA A', 'GSA D', 'GSA E',]
#                'Shannon']

opt_models = ['ohara-2011', 'model-list']
opt_measures = ['LSA-A', 'LSA-D', 'LSA-E', 'GSA-A', 'GSA-D', 'GSA-E',]
#                'Shannon']

savedir = './fig'
if not os.path.isdir(savedir):
    os.makedirs(savedir)


# Create model
model = method.model.VCModel(
    '../mmt/ohara-2011.mmt',
    transform=None,
    dt=dt,
    n_steps=n_steps,
)

# Model parameters
parameters = np.ones(model.n_parameters())


# Get optimal protocols and plot
fig, axes = plt.subplots(len(opt_measures), len(opt_models), sharey=True,
                         figsize=(7, 8))

for i, opt_model in enumerate(opt_models):
    for j, opt_measure in enumerate(opt_measures):
        # Add titles
        if j == 0:
            axes[j, i].set_title(model_side[i] + '\n')

        # Remove x ticks
        #if j != len(opt_measures) - 1:
        axes[j, i].tick_params(axis='x', labelbottom=False)

        # Load protocol
        opt_file = '../design/out/' + opt_measure + '-vc-' + opt_model \
                + '/opt-prt-run0-rank0.txt'
        try:
            all_p = np.loadtxt(opt_file)
        except OSError:
            continue

        # Reshape it to [step_1_voltage, step_1_duration, ...]
        all_p = all_p.flatten().round()
        times = np.arange(0, np.sum(all_p[1::2]), dt)  # ms

        # Set time ticks
        axes[j, i].set_xticks(np.arange(times[0], times[-1], 5000))  # 5 s

        # Update protocol
        model.set_voltage_protocol(all_p)

        # Plot
        color = '#7f7f7f'
        axes[j, i].plot(times, model.voltage(times), c=color)
        axes[j, i].tick_params(axis='y', labelcolor=color)

        # Twinx: current
        ax2 = axes[j, i].twinx()

        color = 'C3'
        ax2.plot(times, model.simulate(parameters, times), c=color)
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.set_ylim((-5, 3))
        if i != len(opt_models) - 1:
            ax2.tick_params(axis='y', labelright=False)
        elif j == 3:
            ax2.set_ylabel('Current\n(A/F)', color='C3')

axes[-1, 1].set_xlabel('\nTime (5000 ms)')
axes[3, 0].set_ylabel('Voltage\n(mV)', color='#7f7f7f')
for j in range(len(opt_measures)):
    axes[j, 0].text(-0.5, 0.5, measure_side[j],
                    transform=axes[j, 0].transAxes, ha='center', va='center',
                    rotation=90)
fig.tight_layout()
fig.savefig('%s/fig1.pdf' % savedir, format='pdf')
plt.close('all')
