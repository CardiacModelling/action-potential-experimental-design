#!/usr/bin/env python3
"""
# Experimental esign for voltage-clamp experiments.
"""
import sys
sys.path.append('..')
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

import psutil
print(psutil.virtual_memory())

import method.model

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
parser.add_argument('-l', '--model_file', help='A mmt model file.')
parser.add_argument('-n', '--n_optim', type=int, default=3,
    help='Number of optimisation repeats.')
parser.add_argument('-r', '--repeat_id', type=int, default=0,
    help='Repeat ID for `--tmp`, splitting a full run into multiple runs.')
parser.add_argument('-p', '--parallel', type=int, default=-1,
    help='Enables/disables parallel evaluation for PINTS.' \
	 + ' Set -1 to use all CPU. Set 0 to disable parallel evaluation.')
parser.add_argument('--rand_x0', action='store_true',
    help='Randomly initialise the optimisation starting point.')
parser.add_argument('--tmp', action='store_true',
    help='Output to tmp folder before postprocessing.')
parser.add_argument('--debug', action='store_true')
args = parser.parse_args()

model_name = os.path.splitext(os.path.basename(args.model_file))[0]

if args.parallel == -1:
    set_parallel = True
else:
    set_parallel = args.parallel

savedir = './out/' + args.design + '-vc-' + model_name
if args.tmp:
    savedir += '-tmp'
    seed_id += args.repeat_id  # Each repeat run has a different seed.
    run_id = str(run_id) + '-re' + str(args.repeat_id)

if not os.path.isdir(savedir):
    os.makedirs(savedir)

np.random.seed(seed_id)
print('Seed ID: ', seed_id)

# Design settings
if 'LSA' in args.design:
    transform = None
    h = 1e-3
    default_param = np.ones(len(method.model.parameters))
elif 'GSA' in args.design or 'Shannon' in args.design:
    transform = np.exp
    n_samples = 512
    logp_lower = [-2] * len(method.model.parameters)  # maybe +/-3
    logp_upper = [2] * len(method.model.parameters)

# Create model
model = method.model.VCModel(
    args.model_file,
    transform=transform,
    dt=dt,
    n_steps=n_steps
)

# Protocol parameter: [step_1_voltage, step_1_duration, step_2..., step3...]
lower = [-120, 50] * n_steps
upper = [60, 2e3] * n_steps
boundaries = pints.RectangularBoundaries(lower, upper)
transformation = pints.RectangularBoundariesTransformation(boundaries)

# Create design
d, c = design_list[args.design]
if 'LSA' in args.design:
    method = pyoed.CentralDifferenceSensitivity
    method_kw = dict(h=h)
    design = d(model, default_param, criterion=c, method=method,
               method_kw=method_kw)
    design.set_maximum_error(1e16)
elif 'GSA' in args.design:
    method = pyoed.SobolFirstOrderSensitivity
    method_kw = dict(n_samples=n_samples)
    b = np.array([logp_lower, logp_upper]).T
    design = d(model, b, criterion=c, method=method, method_kw=method_kw)
    design.set_n_batches(int(n_samples / 2**6))
    design.set_maximum_error(1e16)
elif 'Shannon' in args.design:
    design = None
    raise NotImplementedError

p_evaluate = np.copy(design._method.ask())

# DEBUG: Test parameter samples with a simple protocol
if args.debug:
    test_prt = [-80, 200, 20, 500, -40, 500, -80, 200]
    test_t = np.arange(0, np.sum(test_prt[1::2]), dt)
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
    f.write('\nrand_x0 = ' + str(args.rand_x0))
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

for i_optim in range(args.n_optim):
    # Get x0
    if args.rand_x0:
        need_x0 = True
        while need_x0:
            x0 = boundaries.sample(1)[0]
            if np.isfinite(design(x0)):
                need_x0 = False
    else:
        x0 = [-80, 200, 20, 500, -40, 500, -80, 500] * (n_steps // 4)
    print('x0: ', x0)

    # Try it with x0
    print('Score at x0:', design(x0))
    for _ in range(3):
        assert(design(x0) == design(x0))
    
    opt = pints.OptimisationController(
            design,
            x0,
            sigma0=[5, 50] * n_steps,  # Don't let it make a big jump first.
            boundaries=boundaries,
            transformation=transformation,
            method=pints.CMAES)
    # if set_parallel > opt.optimiser().suggested_population_size():
    #     print('Suggested population size was:',
    #           opt.optimiser().suggested_population_size())
    #     opt.optimiser().set_population_size(set_parallel)
    opt.optimiser().set_population_size(30)
    opt.set_max_iterations(1000)
    opt.set_max_unchanged_iterations(iterations=100, threshold=1e-3)
    opt.set_parallel(set_parallel)

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

    if args.tmp:
        fn = '%s/%s-run%s-%s.txt' % (savedir, prefix, run_id, i_optim)
        method.utils.save_protocol(fn, p)
        gn = '%s/%s-score-run%s-%s.txt' % (savedir, prefix, run_id, i_optim)
        np.savetxt(gn, np.asarray(s).reshape(-1))

    del(opt)
    gc.collect()

if not args.tmp:
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
        fn = '%s/%s-run%s-rank%s.txt' % (savedir, prefix, run_id, i)
        method.utils.save_protocol(fn, obtained_parameters[i])
        gn = '%s/%s-score-run%s-rank%s.txt' % (savedir, prefix, run_id, i)
        np.savetxt(gn, np.asarray(obtained_scores[i]).reshape(-1))

print('Done')
