#!/usr/bin/env python3
"""
# Experimental esign for current-clamp experiments.
"""
from __future__ import print_function
import sys
sys.path.append('../method')
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pints
import pyoed

import model as m

prefix = 'opt-prt'
n_steps = 10  # number of steps of the protocol
dt = 0.1  # ms
seed_id = 101  # random seed
run_id = 0

design_list = {
    'LSA-A':(pyoed.LocalSensitivityDesignMeasure, pyoed.A_criterion),
    'LSA-D':(pyoed.LocalSensitivityDesignMeasure, pyoed.D_criterion),
    'LSA-E':(pyoed.LocalSensitivityDesignMeasure, pyoed.Estar_criterion),
    'GSA-A':(pyoed.GlobalSensitivityDesignMeasure, pyoed.A_criterion),
    'GSA-D':(pyoed.GlobalSensitivityDesignMeasure, pyoed.D_criterion),
    'GSA-E':(pyoed.GlobalSensitivityDesignMeasure, pyoed.Estar_criterion),
    'Shannon':None, # TODO
}

parser = argparse.ArgumentParser('OED for CC experiments.')
parser.add_argument('-d', '--design', type=str,
    choices=design_list.keys(), help='Design for OED.')
parser.add_argument('-l', '--model_file', help='A mmt model file.')
parser.add_argument('-n', '--n_optim', type=int, default=3,
    help='Number of optimisation repeats.')
parser.add_argument('--debug', action='store_true')
args = parser.parse_args()

model_name = os.path.splitext(os.path.basename(args.model_file))[0]

savedir = './out/' + args.design + '-cc-' + model_name
if not os.path.isdir(savedir):
    os.makedirs(savedir)

np.random.seed(seed_id)
print('Seed ID: ', seed_id)

# Design settings
if 'LSA' in args.design:
    transform = None
    h = 1e-3
    default_param = np.ones(len(m.parameters))
elif 'GSA' in args.design or 'Shannon' in args.design:
    transform = np.exp
    n_samples = 1024
    logp_lower = [-2] * len(m.parameters)  # maybe +/-3
    logp_upper = [2] * len(m.parameters)

# Create model
model = m.CCModel(args.model_file, transform=transform, dt=dt, n_steps=n_steps)

# Protocol parameter: [stim_1_amp, holding_1_duration, stim_2..., stim_3...]
lower = [0, 50] * n_steps
upper = [2, 2e3] * n_steps
boundaries = pints.RectangularBoundaries(lower, upper)

# Create design
d, c = design_list[args.design]
if 'LSA' in args.design:
    method = pyoed.CentralDifferenceSensitivity
    method_kw = dict(h=h)
    design = d(model, default_param, criterion=c, method=method,
               method_kw=method_kw)
elif 'GSA' in args.design:
    method = pyoed.SobolFirstOrderSensitivity
    method_kw = dict(n_samples=n_samples)
    b = np.array([logp_lower, logp_upper]).T
    design = d(model, b, criterion=c, method=method, method_kw=method_kw)
elif 'Shannon' in args.design:
    design = None
    raise NotImplementedError
    
p_evaluate = np.copy(design._method.ask())

# DEBUG: Test parameter samples with a simple protocol
if args.debug:
    test_prt = [1, 1000, 2, 800, 1, 500, 0.5, 500]
    test_t = np.arange(0, np.sum(test_prt[1::2]) + 4, dt)
    model.design(test_prt)
    for p in p_evaluate:
        plt.plot(test_t, model.simulate(p))
    plt.xlabel('Times (ms)')
    plt.ylabel('Current (pA)')
    plt.savefig('%s/run%s-test' % (savedir, run_id))
    plt.close()

# DEBUG: Time the evaluation of the design
if args.debug:
    # import timeit
    # print('Single score evaluation time: %s s' \
    #     % (timeit.timeit(lambda: design(x0), number=10) / 10.))
    import cProfile
    import pstats
    fname = savedir + '/design-profile.prof'
    print('Benchmarking design.__call__')
    # Testing protocol parameters
    x0 = [1, 1000, 2, 800] * 5
    cProfile.run('x = design(x0)', fname)
    p = pstats.Stats(fname)
    p.strip_dirs()
    p.sort_stats('time').print_stats(15)  # or 'cumulative'
    print('Calculated design score: ', x)

# Control fitting seed
fit_seed = np.random.randint(0, 2**30)
print('Fit seed: ', fit_seed)
np.random.seed(fit_seed)

# Save parameter for estimating derivative
p_evaluate_file = '%s-run%s-parameter_evaluate.txt' % (prefix, run_id)
np.savetxt(savedir + '/' + p_evaluate_file, p_evaluate)

# Log inputs
with open('%s/%s-run%s.out' % (savedir, prefix, run_id), 'w') as f:
    f.write('design = ' + str(args.design))
    f.write('\nrun_id = ' + str(run_id))
    f.write('\nModel file:')
    f.write('\n    ' + args.model_file)
    f.write('\ntransform = ' + str(transform))
    if 'LSA' in args.design:
        f.write('\ndefault_param = (' \
                + (',').join([str(p) for p in default_param]) + ')')
        f.write('\nh = ' + str(h))
    elif 'GSA' in args.design or 'Shannon' in args.design:
        f.write('\nModel lower bound:\n    ')
        for l in logp_lower:
            f.write(str(l) + ', ')
        f.write('\nModel upper bound:\n    ')
        for u in logp_upper:
            f.write(str(u) + ', ')
        f.write('\nn_samples = ' + str(n_samples))
    f.write('\np_evaluate_file = "' + p_evaluate_file + '"')
    f.write('\nn_steps = ' + str(n_steps))
    f.write('\ndt = ' + str(dt))
    f.write('\nseed_id = ' + str(seed_id))
    f.write('\nfit_seed = ' + str(fit_seed))
    f.write('\nProtocol lower bound:\n    ')
    for l in lower:
        f.write(str(l) + ', ')
    f.write('\nProtocol upper bound:\n    ')
    for u in upper:
        f.write(str(u) + ', ')

# Optimise design
params, scores = [], []

for _ in range(args.n_optim):
    # Get x0
    need_x0 = True
    while need_x0:
        x0 = boundaries.sample(1)[0]
        if np.isfinite(design(x0)):
            need_x0 = False
    print('x0: ', x0)

    # Try it with x0
    print('Score at x0:', design(x0))
    for _ in range(3):
        assert(design(x0) == design(x0))
    
    opt = pints.OptimisationController(
            design,
            x0,
            boundaries=boundaries,
            method=pints.CMAES)
    opt.set_max_iterations(None)
    opt.set_max_unchanged_iterations(iterations=100, threshold=1e-3)
    opt.set_parallel(True)

    # Run optimisation
    try:
        # Tell numpy not to issue warnings
        with np.errstate(all='ignore'):
            p, s = opt.run()
            params.append(p)
            scores.append(s)
            print('Found solution:' )
            print('Voltage (mV)\tDuration (ms)' )
            for i in range(n_steps):
                print(pints.strfloat(p[2 * i]) + '\t' +
                        pints.strfloat(p[2 * i + 1]))
    except ValueError:
        import traceback
        traceback.print_exc()
        raise RuntimeError('Not here...')

# Order from best to worst
order = np.argsort(scores)  # (use [::-1] for LL)
scores = np.asarray(scores)[order]
params = np.asarray(params)[order]

# Show results
bestn = min(5, args.n_optim)
print('Best %d scores:' % bestn)
for i in range(bestn):
    print(scores[i])
print('Mean & std of logposterior:')
print(np.mean(scores))
print(np.std(scores))
print('Worst logposterior:')
print(scores[-1])

#
# Store bestn results
#
obtained_scores = scores[:3]
obtained_parameters = params[:3]
for i in range(bestn):
    p = obtained_parameters[i]
    fn = '%s/%s-run%s-rank%s.txt' % (savedir, prefix, run_id, i)
    with open(fn, 'w') as f:
        f.write('# Voltage [mV]\tDuration [ms]\n')
        for i in range(len(p) // 2):
            f.write(pints.strfloat(p[2 * i]) \
                    + '\t' \
                    + pints.strfloat(p[2 * i + 1]) \
                    + '\n' \
                    )

print('Done')
