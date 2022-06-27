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
dt = 0.1  # ms

model_side = ['Single model', 'Averaged model']
measure_side = ['LSA A', 'LSA D', 'LSA E',]
#                'GSA A', 'GSA D', 'GSA E',]
#                'Shannon']

opt_models = ['ohara-2011', 'model-list']
opt_measures = ['LSA-A', 'LSA-D', 'LSA-E',]
#                'GSA-A', 'GSA-D', 'GSA-E',]
#                'Shannon']

savedir = './fig'
if not os.path.isdir(savedir):
    os.makedirs(savedir)


# Create model
model = method.model.CCModel(
    '../mmt/ohara-2011.mmt',
    transform=None,
    dt=dt,
    n_steps=n_steps,
    #max_evaluation_time=10,
)

# Model parameters
parameters = np.ones(model.n_parameters())


# Get optimal protocols and plot
fig, axes = plt.subplots(len(opt_measures) * len(opt_models) + 1, 1, sharey=True,
                         figsize=(7, 4.5))

for i, opt_model in enumerate(opt_models):
    for j, opt_measure in enumerate(opt_measures):
        k = i * (len(opt_measures) + 1) + j
        # Add titles
        if j == 0:
            axes[k].set_title(model_side[i] + '\n')

        # Remove x ticks
        #if j != len(opt_measures) - 1:
        axes[k].tick_params(axis='x', labelbottom=False)

        # Load protocol
        opt_file = '../design/out/' + opt_measure + '-cc-' + opt_model \
                + '/opt-prt-run0-rank0.txt'
        try:
            all_p = np.loadtxt(opt_file)
        except OSError:
            continue

        # Reshape and round it
        all_p = all_p.flatten().round()

        # Update protocol
        model.design(all_p)
        times = model.times()  # ms

        # Set time ticks
        axes[k].set_xticks(np.arange(times[0], times[-1], 5000))  # 5 s

        # Plot
        color = '#7f7f7f'
        stim = np.append(50, all_p)
        stim += model._stim_dur
        stim = np.cumsum(stim)
        stim_c = np.zeros(times.shape)
        #for s in stim:
        #    #print(np.max(times), np.min(times), len(times), s, model._stim_dur)
        #    stim_c[((times > (s - model._stim_dur)) & (times < s))] = 1
        #axes[k].plot(times, stim_c, '-', c=color)
        for s in stim:
            axes[k].plot([s - model._stim_dur] * 2, [0.1, 0.7], c=color)
        axes[k].tick_params(axis='y', labelcolor=color)
        
        axes[k].set_ylim((-2.5, 1.1))
        axes[k].set_yticks([])

        # Twinx: current
        ax2 = axes[k].twinx()
        s = model.simulate(parameters, times)
        #print(all_p)
        #print(s)

        color = 'C3'
        ax2.plot(times, s, c=color)
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.set_ylim((-105, 165))
        #ax2.set_ylim((-105, 455))
        #if i != len(opt_models) - 1:
        #    ax2.tick_params(axis='y', labelright=False)
        if k in [1, 5]:
            ax2.set_ylabel('Voltage (mV)', color='C3')
            axes[k].set_ylabel('Pacing', color='#7f7f7f')


        axes[k].text(-0.1, 0.5, measure_side[j],
                    transform=axes[k].transAxes, ha='center', va='center',
                    rotation=90)
axes[3].axis('off')
axes[-1].set_xlabel('\nTime (each tick = 5s)')
#fig.tight_layout()
fig.savefig('%s/fig1-cc-v2.pdf' % savedir, format='pdf')
plt.close('all')
