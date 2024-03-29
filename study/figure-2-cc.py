#!/usr/bin/env python3
from __future__ import print_function
import sys
sys.path.append('..')
import os
import numpy as np
import matplotlib.pyplot as plt

import method.heatmap as heatmap  # heap map helper function

import seaborn as sns
cmap = sns.color_palette('crest', as_cmap=True)

#
# Plot all scores for cross criteria for each protocol.
#


# Settings
model_side = ['TNNP', 'Fink', 'OHara', 'Dutta', 'Tomek', 'Averaged']
measure_side = ['LSA A', 'LSA D', 'LSA E']
#, 'GSA A', 'GSA D', 'GSA E',]
#                'Shannon']

opt_models = ['ohara-2011', 'model-list']
opt_measures = ['LSA-A', 'LSA-D', 'LSA-E']
#, 'GSA-A', 'GSA-D', 'GSA-E',]
#                'Shannon']

opt_model_labels = ['Single model', 'Average model']


savedir = './fig/fig2-cc'
if not os.path.isdir(savedir):
    os.makedirs(savedir)


# Get the scores
score_matrix = []
for opt_model in opt_models:
    score_list = []
    for opt_measure in opt_measures:
        # Load score
        score_file = './cross-criteria-evaluate-cc/score-%s-%s.txt' % \
                (opt_model, opt_measure)
        try:
            score = np.loadtxt(score_file)
            score_list.append(score)
        except: # OSError
            score = np.full((len(model_side), len(measure_side)), np.nan)
            score_list.append(score)
    score_matrix.append(score_list)
score_matrix = np.array(score_matrix)

# Benchmark
benchmark = []
benchmark_name = ['biomarkers-vm', '1hz', 'groenendaal-2015-cc']
for b in benchmark_name:
    score_file = './cross-criteria-evaluate-cc/score-%s.txt' % b
    try:
        score = np.loadtxt(score_file)
    except:
        score = np.full((len(model_side), len(measure_side)), np.nan)
    benchmark.append(score)
benchmark = np.array(benchmark)

# Rank the score across all optimal protocols
score_rank = np.full(score_matrix.shape, np.nan)
bench_rank = np.full(benchmark.shape, np.nan)
for i in range(len(model_side)):
    for j in range(len(measure_side)):
        to_be_ranked = score_matrix[:, :, i, j]
        # Use percentage to do scoring instead of the ranking
        # This should highlight the outlier if any!
        finite_tbr = to_be_ranked[np.isfinite(to_be_ranked)]
        mi, ma = np.nanmin(finite_tbr), np.nanmax(finite_tbr)
        ranked = 100. * (to_be_ranked - mi) / (ma - mi)  # lower the better
        score_rank[:, :, i, j] = ranked[:]

        # Benchmark: use the same mi, ma, but not contribute to mi, ma
        ranked = 100. * (benchmark[:, i, j] - mi) / (ma - mi)
        bench_rank[:, i, j] = ranked[:]


# Plot score matrix for each optimal protocol
import matplotlib.patches as patches
score_rank[~np.isfinite(score_rank)] = 200.  # TODO think about how to handle
for i, opt_model in enumerate(opt_models):
    for j, opt_measure in enumerate(opt_measures):
        fig, axes = plt.subplots()
        im, cbar = heatmap.heatmap(score_rank[i, j], model_side, measure_side,
                ax=axes, cmap=cmap, cbar_kw=dict(extend='max'),
                cbarlabel='Normalised score (%)')
        clim = (0, 100.)
        thres = (clim[1] - clim[0]) * 0.2 + clim[0]
        im.set_clim(clim)
        texts = heatmap.annotate_heatmap(im, valfmt='{x:.1f}', threshold=thres)
        y_opt = (1.5 + i * 3) + 0.05
        x_opt = (j - 0.5) + 0.05
        axes.add_patch(patches.Rectangle(
            (x_opt, y_opt),
            1. - 0.1, 1. - 0.1,
            edgecolor='red', fill=False, lw=2
        ))
        fig.tight_layout()
        fig.savefig('%s/opt-protocol-%s-%s.pdf' % (savedir, opt_model,
            opt_measure), format='pdf', bbox_inches='tight')
        plt.close(fig)

for i, bench in enumerate(benchmark_name):
    fig, axes = plt.subplots()
    im, cbar = heatmap.heatmap(bench_rank[i], model_side, measure_side,
            ax=axes, cmap=cmap, cbar_kw=dict(extend='max'),
            cbarlabel='Normalised score (%)')
    clim = (0, 100.)
    thres = (clim[1] - clim[0]) * 0.2 + clim[0]
    im.set_clim(clim)
    texts = heatmap.annotate_heatmap(im, valfmt='{x:.1f}', threshold=thres,
            fontsize=8)
    fig.tight_layout()
    fig.savefig('%s/opt-protocol-%s.pdf' % (savedir, bench), format='pdf',
            bbox_inches='tight')
    plt.close(fig)


# Plot averaged score matrix
averaged_score_rank = np.mean(score_rank, axis=(2, 3))  # over measures
averaged_bench_rank = np.mean(bench_rank, axis=(1, 2)).reshape(1, -1)
clim = (np.min(averaged_score_rank), np.max(averaged_score_rank))
thres = (clim[1] - clim[0]) * 0.2 + clim[0]

fig, axes = plt.subplots(2, 1, gridspec_kw={'height_ratios': [6, 3]},
                               figsize=(6, 4))
ax1, ax2 = axes

im1, _ = heatmap.heatmap(averaged_score_rank, opt_model_labels, measure_side,
                           ax=ax1, cmap=cmap,
                           #cbarlabel='Averaged ranking score (%)')
                           cbarlabel=None)
_ = heatmap.annotate_heatmap(im1, valfmt='{x:.1f}', threshold=thres)
im1.set_clim(clim)

names = ['O', 'P', 'Q', 'R', 'S', 'T']
for i, name in enumerate(names):
    xi, yi = i % 3, int(i / 3)
    x = xi * (1 / 3.) + (1 / 3.) * 0.05
    y = (2 - yi) * (1 / 2.) - (1 / 2.) * 0.05
    color = 'white' if i in [1, 2, 4, 5] else 'black'
    ax1.text(x, y, name, transform=ax1.transAxes,
             ha='left', va='top', weight='bold', color=color)

im2, _ = heatmap.heatmap(averaged_bench_rank,
                         ['Benchmark'],
                         ['Biomarkers', '1 Hz', 'Groenendaal et al.'],
                         ax=ax2, cmap=cmap, cbarlabel=None,
                         rotation=-15)
_ = heatmap.annotate_heatmap(im2, valfmt='{x:.1f}', threshold=thres)
im2.set_clim(clim)

names = ['U', 'V', 'W']
for i, name in enumerate(names):
    xi, yi = i % 3, int(i / 3)
    x = xi * (1 / 3.) + (1 / 3.) * 0.05
    y = (1 - yi) * (1 / 1.) - (1 / 1.) * 0.05
    color = 'white'
    ax2.text(x, y, name, transform=ax2.transAxes,
             ha='left', va='top', weight='bold', color=color)

fig.tight_layout()

#fig.subplots_adjust(right=0.9)
#cbar_ax = fig.add_axes([0.925, 0.15, 0.05, 0.8])
cbar = fig.colorbar(im1, ax=axes.ravel().tolist())
cbar.ax.set_ylabel('Averaged cross-measures (%)', rotation=-90, va="bottom")

fig.savefig('%s/opt-protocol-averaged.pdf' % (savedir), format='pdf',
        bbox_inches='tight')
plt.close('all')
