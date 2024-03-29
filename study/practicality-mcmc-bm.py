#!/usr/bin/env python3
import sys
sys.path.append('..')
import os
import numpy as np
import matplotlib.pyplot as plt
import pints
import pints.io
import pints.plot

import method.model

#
# A practical assessment of the optimal protocol. We try to fit the models to
# the optimal protocol generated synthetic data (either with the same model or
# with a different model), and check how good the inferred models/fits are.
#


try:
    input_file = sys.argv[1]
except:
    print('Usage: python %s [str:input_file]' % os.path.basename(__file__))
    sys.exit()


# Get all input variables
import importlib
sys.path.append(os.path.dirname(input_file))
base = os.path.basename(input_file)
info = importlib.import_module(os.path.splitext(base)[0])

print('Run ID: ', info.run_id)


# Settings
dt = 0.1  # ms
seed_id = 101  # random seed
np.random.seed(seed_id)
fit_seed = np.random.randint(0, 2**30)
np.random.seed(fit_seed)
n_steps = 0

model_side = ['TNNP', 'Fink', 'Grandi', 'O\'Hara', 'CiPA', 'Tomek']
measure_side = ['LSA A', 'LSA D', 'LSA E', 'GSA A', 'GSA D', 'GSA E',]
#                'Shannon']

opt_models = ['ohara-2011', 'model-list']
opt_measures = ['LSA-A', 'LSA-D', 'LSA-E', 'GSA-A', 'GSA-D', 'GSA-E',]
#                'Shannon']

savedir = './practicality-mcmc-bm'
if not os.path.isdir(savedir):
    os.makedirs(savedir)


# Load protocol
assert(info.protocol_file == 'biomarkers')


# Create true model and synthetic data
model_true = method.model.CCBiomarkerModel(info.true_model_file, transform=None, dt=dt)
times = model_true.times()
parameters_true = np.ones(model_true.n_parameters())
raw_data = model_true._simulate(parameters_true)
raw_data += np.random.normal(0, info.noise_sigma, size=raw_data.shape)
model_true._biomarkers.set_data(model_true._pacing, times, raw_data)
data = model_true.extract()
if '--debug' in sys.argv:
    print(model_true._biomarkers.list_of_biomarkers)
    print(data)
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 1, figsize=(6,4), sharex=True)
    axes.plot(times, raw_data)
    axes.set_ylabel('Voltage (mV)')
    axes.set_xlabel('Time (ms)')
    plt.show()
    sys.exit()
del(model_true)


# Estimate noise
noise_sigma = np.std(raw_data[:500])
print('Estimated noise sigma: ', noise_sigma)


# Create fitting model
model_fit = method.model.CCBiomarkerModel(info.fit_model_file, transform=None, dt=dt)
if '--debug2' in sys.argv:
    print(model_fit._biomarkers.list_of_biomarkers)
    print(model_fit.simulate(parameters_true))
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 1, figsize=(6,4), sharex=True)
    axes.plot(times, raw_data)
    axes.plot(times, model_fit._simulate(parameters_true))
    axes.set_ylabel('Voltage (mV)')
    axes.set_xlabel('Time (ms)')
    plt.show()
    sys.exit()

# Create likelihood
model_fit.set_denominators(np.array(data))
problem = pints.SingleOutputProblem(model_fit, np.arange(len(model_fit._biomarkers.list_of_biomarkers)), np.ones(len(data)))
# log_likelihood = pints.GaussianLogLikelihood(problem)
log_likelihood = pints.GaussianKnownSigmaLogLikelihood(problem, noise_sigma / np.std(raw_data))
log_prior = pints.ComposedLogPrior(
    *([pints.LogNormalLogPrior(0, 0.15)] * (problem.n_parameters()))# + 1))
)
log_posterior = pints.LogPosterior(log_likelihood, log_prior)

transformation = pints.LogTransformation(problem.n_parameters())# + 1)

# Check log_posterior is not throwing error and deterministic
for _ in range(3):
    assert(log_posterior(np.ones(model_fit.n_parameters())) ==\
            log_posterior(np.ones(model_fit.n_parameters())))# + 1)))

# Get samples with finite posteriors from priors
mcmc_init = []
for i in range(4):
    p = np.nan
    while not np.isfinite(p):
        x = log_prior.sample(n=1)[0]
        p = log_posterior(x)
        print(x, p)
    mcmc_init.append(x)
mcmc_init = np.array(mcmc_init)

# Run
saveas = savedir + '/run_%s' % info.run_id
mcmc = pints.MCMCController(log_posterior, 4, mcmc_init,
                            transformation=transformation,
                            method=pints.PopulationMCMC,)
n_iter = 40000
mcmc.set_max_iterations(n_iter)
mcmc.set_initial_phase_iterations(200)
mcmc.set_parallel(True)
mcmc.set_chain_filename('%s-chain.csv' % saveas)
mcmc.set_log_pdf_filename('%s-log-posterior.csv' % saveas)
mcmc.set_log_interval(iters=200, warm_up=5)
chains = mcmc.run()

# Save
pints.io.save_samples('%s-chain.csv' % saveas, *chains)


# Simple plotting of results

# burn in and thinning
chains_final = chains[:, int(0.8 * n_iter)::2, :]

# Plot
pints.plot.pairwise(chains_final[0], kde=False, ref_parameters=parameters_true)#np.append(parameters_true, 1))
plt.savefig('%s-fig1.png' % saveas)
plt.close('all')

pints.plot.trace(chains_final, ref_parameters=parameters_true)#np.append(parameters_true, 1))
plt.savefig('%s-fig2.png' % saveas)
plt.close('all')

#pints.plot.series(chains_final[0], problem)
#plt.savefig('%s-fig3.png' % saveas)
#plt.close('all')

# Check convergence using rhat criterion
print('R-hat:')
r_hats = pints.rhat(chains_final)
print(r_hats)
np.savetxt('%s-rhat.txt' % saveas, r_hats)
