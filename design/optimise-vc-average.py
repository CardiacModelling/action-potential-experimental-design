#!/usr/bin/env python3
"""
# Experimental esign for voltage-clamp experiments with model averaging.
"""
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
import gc
gc.enable()
# gc.set_debug(gc.DEBUG_LEAK)

import model as m

prefix = 'opt-prt'
n_steps = 20  # number of steps of the protocol
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

parser = argparse.ArgumentParser('OED for VC experiments.')
parser.add_argument('-d', '--design', type=str,
    choices=design_list.keys(), help='Design for OED.')
parser.add_argument('-l', '--model_file_list',
    help='a txt file containing a list of mmt file names.')
parser.add_argument('-n', '--n_optim', type=int, default=3,
    help='Number of optimisation repeats.')
parser.add_argument('--debug', action='store_true')
args = parser.parse_args()

file_name = os.path.splitext(os.path.basename(args.model_file_list))[0]

with open(args.model_file_list, 'r') as f:
    ls = f.readlines()
args.model_file_list = [l.strip() for l in ls]

savedir = './out/' + args.design + '-vc-' + file_name
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

# Create models
model_list = []
for model_file in args.model_file_list:
    model_list.append(
        m.VCModel(model_file, transform=transform, dt=dt, n_steps=n_steps))

# Protocol parameter: [step_1_voltage, step_1_duration, step_2..., step3...]
lower = [-120, 50] * n_steps
upper = [60, 2e3] * n_steps
boundaries = pints.RectangularBoundaries(lower, upper)

# Create design
d, c = design_list[args.design]
design_list = []
if 'LSA' in args.design:
    method = pyoed.CentralDifferenceSensitivity
    method_kw = dict(h=h)
    for model in model_list:
        design_list.append(
            d(model, default_param, criterion=c, method=method,
              method_kw=method_kw))
elif 'GSA' in args.design:
    method = pyoed.SobolFirstOrderSensitivity
    method_kw = dict(n_samples=n_samples)
    b = np.array([logp_lower, logp_upper]).T
    for model in model_list:
        design = d(model, b, criterion=c, method=method, method_kw=method_kw)
        design.set_n_batches(int(n_samples / 2**8))
        design_list.append(design)
elif 'Shannon' in args.design:
    design = None
    raise NotImplementedError
design = pyoed.CombineDesignMeasure(design_list, aggregate=np.mean)

p_evaluate = np.copy(design._measures[0]._method.ask())

# DEBUG: Test parameter samples with a simple protocol
if args.debug:
    test_prt = [-80, 200, 20, 500, -40, 500, -80, 200]
    test_t = np.arange(0, np.sum(test_prt[1::2]), dt)
    for model in model_list:
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
    x0 = [-80, 200, 20, 500, -40, 500, -80, 500] * (n_steps // 4)
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
    f.write('\nModel files:')
    for model_file in args.model_file_list:
        f.write('\n    ' + model_file)
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
    opt.optimiser().set_population_size(30)
    opt.set_max_iterations(1000)
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
    del(opt)
    gc.collect()

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
obtained_scores = scores[:bestn]
obtained_parameters = params[:bestn]
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
