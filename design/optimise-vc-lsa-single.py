#!/usr/bin/env python3
"""
# LSA-based design for voltage-clamp experiments.
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
n_steps = 20  # number of steps of the protocol
dt = 0.1  # ms
seed_id = 101  # random seed
run_id = 0

criterion_list = dict(
    A=pyoed.A_criterion,
    D=pyoed.D_criterion,
    E=pyoed.Estar_criterion,
)

parser = argparse.ArgumentParser('LSA for VC experiments.')
parser.add_argument('-c', '--criterion', type=str,
    choices=criterion_list.keys(), help='Criterion for LSA design.')
parser.add_argument('-l', '--model_file', help='A mmt model file.')
parser.add_argument('--debug', action='store_true')
args = parser.parse_args()

model_name = os.path.splitext(os.path.basename(args.model_file))[0]

savedir = './out/lsa-' + args.criterion + '-vc-' + model_name
if not os.path.isdir(savedir):
    os.makedirs(savedir)

np.random.seed(seed_id)
print('Seed ID: ', seed_id)

model = m.VCModel(args.model_file, transform=None, dt=dt, n_steps=n_steps)

# protocol parameter: [step_1_voltage, step_1_duration, step_2..., step3...]
lower = [-120, 50] * n_steps
upper = [60, 2e3] * n_steps
boundaries = pints.RectangularBoundaries(lower, upper)

# Create design
default_param = np.ones(model.n_parameters())
h = 1e-3
design = pyoed.LocalSensitivityDesignMeasure(model, default_param,
    criterion=criterion_list[args.criterion])
p_evaluate = np.copy(design._method.ask())

# DEBUG: Test parameter samples with a simple protocol
if args.debug:
    test_prt = [-80, 200, 20, 500, -40, 500, -80, 200]
    test_t = np.arange(0, np.sum(test_prt[1::2]), dt)
    model.set_voltage_protocol(test_prt)
    for p in p_evaluate:
        plt.plot(test_t, model.simulate(p, times=test_t))
    plt.xlabel('Times (ms)')
    plt.ylabel('Current (pA)')
    plt.savefig('%s/%s-run%s-test' % (savedir, prefix, run_id))
    plt.close()

# Simple initial guess of the protocol
x0 = [-80, 200, 20, 500, -40, 500, -80, 500] * (n_steps // 4)

# DEBUG: Time the evaluation of the design
if args.debug:
    #import timeit
    #print('Single score evaluation time: %s s' \
    #    % (timeit.timeit(lambda: design(x0), number=10) / 10.))
    import cProfile
    import pstats
    fname = 'design-lsa-profile.prof'
    print('Benchmarking design.__call__')
    cProfile.run('x = design(x0)', fname)
    p = pstats.Stats(fname)
    p.strip_dirs()
    p.sort_stats('time').print_stats(15)
    #p.sort_stats('cumulative').print_stats(15)
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
    f.write('run_id = ' + str(run_id))
    f.write('\nModel file:')
    f.write('\n    ' + args.model_file)
    f.write('\ndefault_param = (' \
            + (',').join([str(p) for p in default_param]) + ')')
    f.write('\nh = ' + str(h))
    f.write('\np_evaluate_file = "' + p_evaluate_file + '"')
    f.write('\nn_steps = ' + str(n_steps))
    f.write('\ndt = ' + str(dt))
    f.write('\ncriterion = ' + str(args.criterion))
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

N = 3

for _ in range(N):
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
bestn = min(3, N)
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
