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
# Compute eigenvalues.
#


# Settings
model_file_list = [
    #'../mmt/tnnp-2004.mmt',
    #'../mmt/fink-2008.mmt',
    #'../mmt/grandi-2010.mmt',
    '../mmt/ohara-2011.mmt',
    #'../mmt/cipa-2017.mmt',
    #'../mmt/tomek-2019.mmt',
]

class LSDesignMeasure(pyoed.LocalSensitivityDesignMeasure):
    def __call__(self, x):
        # Update design
        self._model.design(x)

        # Get points to calculate sensitivity
        ps = self._method.ask()
        fs = []
        for p in ps:
            fs.append(self._model.simulate(p))
        fs = np.asarray(fs)
        if not np.isfinite(fs).all():
            return self._max_error

        # Calculate sensitivity
        if self._n_batches is None:
            s = self._method.tell(fs)
        else:
            n_outputs = len(fs[0])
            s = np.zeros((self._n_mp, n_outputs))
            s[:] = np.nan
            n_each = n_outputs // self._n_batches
            start = 0
            if n_each > 0:
                for i in range(self._n_batches - 1):
                    end = (i + 1) * n_each
                    s[:, start:end] = self._method.tell(fs[:, start:end])
                    start = end
            s[:, start:] = self._method.tell(fs[:, start:])
            s = s.T  # rows being model outputs; cols being model parameters

        # Calculate measure
        return self._criterion(s)

class GSDesignMeasure(pyoed.GlobalSensitivityDesignMeasure):
    def __call__(self, x):
        # Update design
        self._model.design(x)

        # Get points to calculate sensitivity
        ps = self._method.ask()
        fs = []
        for p in ps:
            fs.append(self._model.simulate(p))
        fs = np.asarray(fs)
        if not np.isfinite(fs).all():
            return self._max_error

        # Calculate sensitivity
        if self._n_batches is None:
            s = self._method.tell(fs)
        else:
            n_outputs = len(fs[0])
            s = np.zeros((self._n_mp, n_outputs))
            s[:] = np.nan
            n_each = n_outputs // self._n_batches
            start = 0
            if n_each > 0:
                for i in range(self._n_batches - 1):
                    end = (i + 1) * n_each
                    s[:, start:end] = self._method.tell(fs[:, start:end])
                    start = end
            s[:, start:] = self._method.tell(fs[:, start:])
            s = s.T  # rows being model outputs; cols being model parameters

        # Calculate measure
        return self._criterion(s)

design_list = OrderedDict(
    LSA=LSDesignMeasure,
    GSA=GSDesignMeasure,
)

n_samples = 512  # number of samples to be compared
n_steps = 20  # number of steps of the protocol
dt = 0.5  # ms
seed_id = 101  # random seed

opt_models = ['ohara-2011', 'model-list']
opt_measures = ['LSA-A', 'LSA-D', 'LSA-E', 'GSA-A', 'GSA-D', 'GSA-E',]
#                'Shannon']

savedir = './eigenvalues-vc'
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

def c(s):
    """
    Compute the eigenvalues for a given sensitivity matrix `s`.
    output = eigenvalues([s^T s]^-1)

    s: the sensitivity matrix
       [s_y1_theta1, s_y1_theta2, ...]
       [s_y2_theta1, s_y2_theta2, ...]
       [...        , ...        , ...].
    """
    try:
        i_matrix = np.linalg.pinv(np.matmul(s.T, s))
        eigvals = np.linalg.eigvals(i_matrix)
        return eigvals
    except np.linalg.LinAlgError:
        return float('nan')

# Create scores shape (model_file_list, design_list)
score_list = []
# Single models
for i_model in range(len(model_list)):
    for n in design_list:
        d = design_list[n]
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

'''
# Average model
score_list = []
for n in design_list:
    average_list = []
    d = design_list[n]
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
'''


# Compute for random protocols
np.random.seed(seed_id)
boundaries = pints.RectangularBoundaries(lower, upper)
for i in range(5):
    x0 = boundaries.sample(1)[0]

    # Loop over LSA and GSA
    score_per_prt = []
    for score in score_list:
        score_per_prt.append(score(x0))

    # Save score matrix
    np.savetxt('%s/eigenvalues-lsa-gsa-rand-%s.txt' % (savedir, i),
            np.array(score_per_prt).T)


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

        # Loop over LSA and GSA
        for score in score_list:
            score_per_prt.append(score(all_p))

        # Save score matrix
        np.savetxt('%s/eigenvalues-lsa-gsa-%s-%s.txt' % (savedir, opt_model, opt_measure),
                np.array(score_per_prt).T)


# Get benchmark protocols
for opt_file_name in ['ch3', 'groenendaal-2015']:
    # Load protocol
    opt_file = './benchmark-protocols/' + opt_file_name + '.txt'
    try:
        all_p = np.loadtxt(opt_file)
    except: # OSError
        continue

    # Reshape it to [step_1_voltage, step_1_duration, ...]
    all_p = all_p.flatten().round()

    # Loop over LSA and GSA
    score_per_prt = []
    for score in score_list:
        score_per_prt.append(score(all_p))

    # Save score matrix
    np.savetxt('%s/eigenvalues-lsa-gsa-%s.txt' % (savedir, opt_file_name),
            np.array(score_per_prt).T)

