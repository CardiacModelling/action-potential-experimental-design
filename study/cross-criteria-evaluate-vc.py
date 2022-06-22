#!/usr/bin/env python3
import sys
sys.path.append('..')
import os
from collections import OrderedDict
import numpy as np
import pints
import pyoed
import method.model

#
# Compute all scores for cross criteria.
#


# Settings
model_file_list = [
    '../mmt/tnnp-2004.mmt',
    '../mmt/fink-2008.mmt',
    #'../mmt/grandi-2010.mmt',
    '../mmt/ohara-2011.mmt',
    '../mmt/cipa-2017.mmt',
    '../mmt/tomek-2019.mmt',
]

design_list = OrderedDict(
    LSA_A=(pyoed.LocalSensitivityDesignMeasure, pyoed.A_criterion),
    LSA_D=(pyoed.LocalSensitivityDesignMeasure, pyoed.D_criterion),
    LSA_E=(pyoed.LocalSensitivityDesignMeasure, pyoed.Estar_criterion),
    GSA_A=(pyoed.GlobalSensitivityDesignMeasure, pyoed.A_criterion),
    GSA_D=(pyoed.GlobalSensitivityDesignMeasure, pyoed.D_criterion),
    GSA_E=(pyoed.GlobalSensitivityDesignMeasure, pyoed.Estar_criterion),
#    Shannon=(ShannonDesignMeasure, [1, 1, 0, 1]),
)

n_samples = 512  # number of samples to be compared
n_steps = 20  # number of steps of the protocol
dt = 0.1  # ms
seed_id = 101  # random seed

opt_models = ['ohara-2011', 'model-list']
opt_measures = ['LSA-A', 'LSA-D', 'LSA-E', 'GSA-A', 'GSA-D', 'GSA-E',]
#                'Shannon']

savedir = './cross-criteria-evaluate-vc'
if not os.path.isdir(savedir):
    os.makedirs(savedir)


# Create models
model_list = []
log_model_list = []
# NOTE: Transform for GSA, not for LSA!
for model_file in model_file_list:
    model_list.append(method.model.VCModel(model_file, transform=None, dt=dt, n_steps=n_steps))
    log_model_list.append(method.model.VCModel(model_file, transform=np.exp, dt=dt, n_steps=n_steps))

# Model parameter bounds
logp_lower = [-2] * len(method.model.parameters)  # maybe +/-3
logp_upper = [2] * len(method.model.parameters)


# protocol parameter: [step_1_voltage, step_1_duration, step_2..., step3...]
lower = [-120, 50] * n_steps
upper = [60, 2e3] * n_steps

# LSA setting
lsa_local_param = np.ones(len(method.model.parameters))
lsa_h = 1e-3
# GSA setting
boundaries = np.array([logp_lower, logp_upper]).T


# Create scores shape (model_file_list, design_list)
score_matrix = []
# Single models
for i_model in range(len(model_list)):
    score_list = []
    for n in design_list:
        d, c = design_list[n]
        if 'LSA' in n:
            model = model_list[i_model]
            sensitivity_method = pyoed.CentralDifferenceSensitivity
            method_kw = dict(h=lsa_h)
            design = d(model, lsa_local_param, criterion=c, method=sensitivity_method,
                       method_kw=method_kw)
        elif 'GSA' in n:
            model = log_model_list[i_model]
            sensitivity_method = pyoed.SobolFirstOrderSensitivity
            method_kw = dict(n_samples=n_samples)
            design = d(model, boundaries, criterion=c, method=sensitivity_method,
                       method_kw=method_kw)
            design.set_n_batches(int(n_samples / 2**8))
        score_list.append(design)
    score_matrix.append(score_list)

# Average model
score_list = []
for n in design_list:
    average_list = []
    d, c = design_list[n]
    if 'LSA' in n:
        sensitivity_method = pyoed.CentralDifferenceSensitivity
        method_kw = dict(h=lsa_h)
        for model in model_list:
            design = d(model, lsa_local_param, criterion=c, method=sensitivity_method,
                       method_kw=method_kw)
            average_list.append(design)
    elif 'GSA' in n:
        sensitivity_method = pyoed.SobolFirstOrderSensitivity
        method_kw = dict(n_samples=n_samples)
        for model in log_model_list:
            design = d(model, boundaries, criterion=c, method=sensitivity_method,
                       method_kw=method_kw)
            design.set_n_batches(int(n_samples / 2**8))
            average_list.append(design)
    score_list.append(pyoed.CombineDesignMeasure(average_list, aggregate=np.mean))
score_matrix.append(score_list)


# Get optimal protocols and evaluate the scores
for opt_model in opt_models:
    for opt_measure in opt_measures:
        # Load protocol
        opt_file = '../design/out/' + opt_measure + '-vc-' + opt_model \
                + '/opt-prt-run0-rank0.txt'
        try:
            all_p = np.loadtxt(opt_file)
        except: # OSError
            continue

        # Reshape it to [step_1_voltage, step_1_duration, ...]
        all_p = all_p.flatten().round()

        # Compute score
        # Score matrix *per optimal protocol*
        score_per_prt = []

        for score_list in score_matrix:
            # Loop over models
            score_per_prt_per_list = []
            for score in score_list:
                score_per_prt_per_list.append(score(all_p))
            score_per_prt.append(score_per_prt_per_list)

        # Save score matrix
        np.savetxt('%s/score-%s-%s.txt' % (savedir, opt_model, opt_measure),
                score_per_prt)
