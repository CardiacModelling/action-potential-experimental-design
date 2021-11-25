#!/usr/bin/env python3
import sys
sys.path.append('..')
import os
import numpy as np
import pints

import model as m
import score_lsa as s_lsa
import score_gsa as s_gsa
import score_shannon as s_sha

#
# Compute all scores for cross criteria.
#


# Settings
model_file_list = [
    '../mmt/tnnp-2004.mmt',
    '../mmt/fink-2008.mmt',
    '../mmt/grandi-2010.mmt',
    '../mmt/ohara-2011.mmt',
    '../mmt/cipa-2017.mmt',
    '../mmt/tomek-2019.mmt',
]

design_list = {
    'LSA-A':(pyoed.LocalSensitivityDesignMeasure, pyoed.A_criterion),
    'LSA-D':(pyoed.LocalSensitivityDesignMeasure, pyoed.D_criterion),
    'LSA-E':(pyoed.LocalSensitivityDesignMeasure, pyoed.Estar_criterion),
    'GSA-A':(pyoed.GlobalSensitivityDesignMeasure, pyoed.A_criterion),
    'GSA-D':(pyoed.GlobalSensitivityDesignMeasure, pyoed.D_criterion),
    'GSA-E':(pyoed.GlobalSensitivityDesignMeasure, pyoed.Estar_criterion),
    'Shannon':(ShannonDesignMeasure, [1, 1, 0, 1]),
}

n_bootstrap_samples = 100  # number of bootstrap model samples
n_samples = 500  # number of samples to be compared
n_steps = 20  # number of steps of the protocol
dt = 0.1  # ms
seed_id = 101  # random seed

model_side = ['TNNP', 'Fink', 'OHara', 'Paci', 'BMA', 'Bootstrap']
measure_side = ['LSA A', 'LSA D', 'LSA E', 'GSA A', 'GSA D', 'GSA E',
        'Shannon']

opt_models = ['', '-bma', '-bs']
opt_measures = ['lsa-A', 'lsa-D', 'lsa-E', 'gsa-sobol-A', 'gsa-sobol-D',
        'gsa-sobol-E', 'shannon']

savedir = './cross-criteria-evaluate'
if not os.path.isdir(savedir):
    os.makedirs(savedir)


# Create models
model_list = []
# Single models
for model_file in model_file_list:
    model_list.append(m.VCModel(model_file, transform=np.exp, dt=dt))
# BMA
model_prior = np.ones(len(model_file_list)) / np.float(len(model_file_list))
model_list.append(m.BMAVCModel(model_file_list, prior=model_prior,
        transform=np.exp, dt=dt))
# Bootstrap
model_list.append(m.BootstrapVCModel(model_file_list,
        n_samples=n_bootstrap_samples, transform=np.exp, dt=dt))

# Model parameter bounds
logp_lower = [-2] * len(m.parameters)  # maybe +/-3
logp_upper = [2] * len(m.parameters)


# protocol parameter: [step_1_voltage, step_1_duration, step_2..., step3...]
lower = [-120, 50] * n_steps
upper = [60, 2e3] * n_steps
boundaries = pints.RectangularBoundaries(lower, upper)


# LSA setting
lsa_local_param = np.ones(model_list[0].n_parameters())
lsa_h = 1e-3
# GSA setting (SALib)
salib_problem = {
    'num_vars': model_list[0].n_parameters(),
    'names': ['p%s' % (i) for i in range(model_list[0].n_parameters())],
    'bounds': np.array([logp_lower, logp_upper]).T,
}
# Shannon measure setting
weight_sha = [
    10.,  # shannon_all_measure
    0.1,  # shannon_segments_measure
    0.,  # parameter_dependency_measure; see [1] why we drop it.
    1.,  # output_uncertainty_measure
]


# Create scores shape (model_side, measure_side)
score_matrix = []
# Single models
for model in model_list[:-2]:
    score_list = []
    # LSA
    for criterion in criterion_list:
        score_list.append(
            s_lsa.LSAScore([model],
                default_param=lsa_local_param,
                h=lsa_h,
                n_steps=n_steps,
                criterion=criterion)
        )
    # GSA
    for criterion in criterion_list:
        score_list.append(
            s_gsa.GSASobolScore([model],
                problem=salib_problem,
                n_steps=n_steps,
                n_samples=n_samples,
                criterion=criterion)
        )
    # Shannon
    score_list.append(
        s_sha.ShannonGSAScore([model],
            problem=salib_problem,
            n_steps=n_steps,
            n_samples=n_samples,
            weight=weight_sha)
    )
    score_matrix.append(score_list)
# BMA
score_list = []
# LSA
for criterion in criterion_list:
    score_list.append(
        s_lsa.LSAScore([model_list[-2]],
            default_param=lsa_local_param,
            h=lsa_h,
            n_steps=n_steps,
            criterion=criterion)
    )
# GSA
for criterion in criterion_list:
    score_list.append(
        s_gsa.GSASobolScore([model_list[-2]],
            problem=salib_problem,
            n_steps=n_steps,
            n_samples=n_samples,
            criterion=criterion)
    )
# Shannon
score_list.append(
    s_sha.ShannonGSAScore([model_list[-2]],
        problem=salib_problem,
        n_steps=n_steps,
        n_samples=n_samples,
        weight=weight_sha)
)
score_matrix.append(score_list)
# Bootstrap
score_list = []
# LSA
for criterion in criterion_list:
    score_list.append(
        s_lsa.LSAScore_Bootstrap(model_list[-1],
            default_param=lsa_local_param,
            h=lsa_h,
            n_steps=n_steps,
            criterion=criterion)
    )
# GSA
for criterion in criterion_list:
    score_list.append(
        s_gsa.GSASobolScore_Bootstrap(model_list[-1],
            problem=salib_problem,
            n_steps=n_steps,
            n_samples=n_samples,
            criterion=criterion)
    )
# Shannon
score_list.append(
    s_sha.ShannonGSAScore_Bootstrap(model_list[-1],
        problem=salib_problem,
        n_steps=n_steps,
        n_samples=n_samples,
        weight=weight_sha)
)
score_matrix.append(score_list)


# Get optimal protocols and evaluate the scores
for opt_model in opt_models:
    for opt_measure in opt_measures:
        # Load protocol
        opt_file = '../design/out-' + opt_measure + '/opt-prt' + opt_model \
                + '-run0-rank0.txt'
        try:
            all_p = np.loadtxt(opt_file)
        except: # OSError
            continue

        # Reshape it to [step_1_voltage, step_1_duration, ...]
        all_p = all_p.flatten()
        times = np.arange(0, np.sum(all_p[1::2]), dt)  # ms

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
        np.savetxt('%s/score%s-%s.txt' % (savedir, opt_model, opt_measure),
                score_per_prt)
