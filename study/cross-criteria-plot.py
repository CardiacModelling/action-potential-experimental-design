#!/usr/bin/env python3
import sys
sys.path.append('..')
import os
import numpy as np
import matplotlib.pyplot as plt
import method.heatmap as heatmap  # heap map helper function

#
# Plot all scores for cross criteria for each protocol.
#


# Settings
model_side = ['TNNP', 'Fink', 'Grandi', 'O\'Hara', 'CiPA', 'Tomek']
measure_side = ['LSA A', 'LSA D', 'LSA E', 'GSA A', 'GSA D', 'GSA E',]
#                'Shannon']

opt_models = ['ohara-2011', 'model-list']
opt_measures = ['LSA-A', 'LSA-D', 'LSA-E', 'GSA-A', 'GSA-D', 'GSA-E',]
#                'Shannon']

opt_model_labels = ['Single model', 'Average model']


savedir = './cross-criteria-plot'
if not os.path.isdir(savedir):
    os.makedirs(savedir)


# Get the scores
score_matrix = []
for opt_model in opt_models:
    score_list = []
    for opt_measure in opt_measures:
        # Load score
        score_file = './cross-criteria-evaluate/score%s-%s.txt' % \
                (opt_model, opt_measure)
        try:
            score = np.loadtxt(score_file)
            score_list.append(score)
        except: # OSError
            score = np.full((len(model_side), len(measure_side)), np.nan)
            score_list.append(score)
    score_matrix.append(score_list)
score_matrix = np.array(score_matrix)


# Rank the score across all optimal protocols
score_rank = np.full(score_matrix.shape, np.nan)
for i in range(len(model_side)):
    for j in range(len(measure_side)):
        to_be_ranked = score_matrix[:, :, i, j]
        # Use percentage to do scoring instead of the ranking
        # This should highlight the outlier if any!
        finite_tbr = to_be_ranked[np.isfinite(to_be_ranked)]
        mi, ma = np.nanmin(finite_tbr), np.nanmax(finite_tbr)
        ranked = 100. * (to_be_ranked - mi) / (ma - mi)  # lower the better
        score_rank[:, :, i, j] = ranked[:]


# Plot score matrix for each optimal protocol
score_rank[~np.isfinite(score_rank)] = 200.  # TODO think about how to handle
for i, opt_model in enumerate(opt_models):
    for j, opt_measure in enumerate(opt_measures):
        fig, axes = plt.subplots()
        im, cbar = heatmap.heatmap(score_rank[i, j], model_side, measure_side,
                ax=axes, cmap='YlGn', cbarlabel='Ranking score (%)')
        clim = (0, 100.)
        thres = (clim[1] - clim[0]) * 0.6
        im.set_clim(clim)
        texts = heatmap.annotate_heatmap(im, valfmt='{x:.1f}', threshold=thres)
        fig.tight_layout()
        fig.savefig('%s/opt-protocol%s-%s' % (savedir, opt_model,
            opt_measure), dpi=200)
        plt.close(fig)


# Plot averaged score matrix
averaged_score_rank = np.mean(score_rank, axis=(2, 3))  # over measures
fig, axes = plt.subplots()
im, cbar = heatmap.heatmap(averaged_score_rank, opt_model_labels, measure_side,
        ax=axes, cmap='YlGn', cbarlabel='Averaged ranking score (%)')
texts = heatmap.annotate_heatmap(im, valfmt='{x:.1f}')
fig.tight_layout()
fig.savefig('%s/opt-protocol-averaged' % (savedir), dpi=200)
plt.close('all')
