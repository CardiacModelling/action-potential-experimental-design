#!/usr/bin/env python3
import numpy as np

# Setting
models = ['tnnp', 'fink', 'grandi', 'ohara', 'cipa', 'tomek']

model_dir = {
    'tnnp': '../mmt/tnnp-2004.mmt',
    'fink': '../mmt/fink-2008.mmt',
    'grandi': '../mmt/grandi-2010.mmt',
    'ohara': '../mmt/ohara-2011.mmt',
    'cipa': '../mmt/cipa-2017.mmt',
    'tomek': '../mmt/tomek-2019.mmt',
}

noise_sigma_vc = {
    'tnnp': 0.5,  # bigger current?
    'fink': 0.15,
    'grandi': 0.15,
    'ohara': 0.15,
    'cipa': 0.15,
    'tomek': 0.15,
}

noise_sigma_cc = {
    'tnnp': 1,
    'fink': 1,
    'grandi': 1,
    'ohara': 1,
    'cipa': 1,
    'tomek': 1,
}

opt_models = ['ohara-2011', 'model-list']
opt_measures = ['LSA-A', 'LSA-D', 'LSA-E', 'GSA-A', 'GSA-D', 'GSA-E',]
#                'Shannon']


def get_prt(opt_measure, opt_model, mode='vc'):
    # Return protocol directory
    i = 0
    prt = '../design/out/' + opt_measure + '-' + mode + '-' + opt_model \
            + '/opt-prt-run0-rank0.txt'
    return prt


def write(f, run_id, true, fit, prt, noise):
    # Write "id_[run_id].py" file
    f.write('run_id = \"' + run_id + '\"\n')
    f.write('true_model_file = \"' + true + '\"\n')
    f.write('fit_model_file = \"' + fit + '\"\n')
    f.write('protocol_file = \"' + prt + '\"\n')
    f.write('noise_sigma = ' + str(noise))


def generate(true, fit, opt_models, opt_measures, int_id, ii='', mode='vc'):
    # Generate "id_[run_id].py" file for all combinations of `opt_models` and
    # `opt_measures` for the given `true` and `fit` models.
    # 
    # Return the updated `int_id`
    ftrue = model_dir[true]
    ffit = model_dir[fit]
    matrix = []
    noise_sigma = noise_sigma_vc if mode == 'vc' else noise_sigma_cc
    for opt_model in opt_models:
        row = []
        for opt_measure in opt_measures:
            run_id = '%03d' % int_id
            row.append(int_id)
            prt = get_prt(opt_measure, opt_model, mode=mode)
            with open('id_%s.py' % run_id, 'w') as f:
                write(f, run_id, ftrue, ffit, prt, noise_sigma[true])
            int_id += 1
        matrix.append(row)
    np.savetxt('true_%s-fit_%s-row_models-col_measures-%s%s.txt' \
            % (true, fit, mode, ii), matrix, fmt='%i')
    return int_id


# Starts
int_id = 1


# Loop through all combinations of true and fit models
for model in models:
    int_id = generate(model, model, opt_models, opt_measures, int_id, mode='vc')


# For benchmark protocols
matrix = []
for true in models:
    m = model_dir[true]
    row = []
    run_id = '%03d' % int_id
    row.append(int_id)
    prt = './benchmark-protocols/ch3.txt'
    with open('id_%s.py' % run_id, 'w') as f:
        write(f, run_id, m, m, prt, noise_sigma_vc[true])
    int_id += 1
    matrix.append(row)
np.savetxt('ch3.txt', matrix, fmt='%i')

matrix = []
for true in models:
    m = model_dir[true]
    row = []
    run_id = '%03d' % int_id
    row.append(int_id)
    prt = './benchmark-protocols/groenendaal-2015.txt'
    with open('id_%s.py' % run_id, 'w') as f:
        write(f, run_id, m, m, prt, noise_sigma_vc[true])
    int_id += 1
    matrix.append(row)
np.savetxt('groenendaal-2015.txt', matrix, fmt='%i')


# Loop through all combinations of true and fit models
for model in models:
    int_id = generate(model, model, opt_models, opt_measures, int_id, mode='cc')

# For benchmark protocols (biomarkers)
matrix = []
for true in models:
    m = model_dir[true]
    row = []
    run_id = '%03d' % int_id
    row.append(int_id)
    prt = 'biomarkers'
    with open('id_%s.py' % run_id, 'w') as f:
        write(f, run_id, m, m, prt, noise_sigma_cc[true])
    int_id += 1
    matrix.append(row)
np.savetxt('biomarkers.txt', matrix, fmt='%i')
