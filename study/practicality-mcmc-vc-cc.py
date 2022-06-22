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
n_steps = 20

model_side = ['TNNP', 'Fink', 'Grandi', 'O\'Hara', 'CiPA', 'Tomek']
measure_side = ['LSA A', 'LSA D', 'LSA E', 'GSA A', 'GSA D', 'GSA E',]
#                'Shannon']

opt_models = ['ohara-2011', 'model-list']
opt_measures = ['LSA-A', 'LSA-D', 'LSA-E', 'GSA-A', 'GSA-D', 'GSA-E',]
#                'Shannon']

savedir = './practicality-mcmc-vc-cc'
if not os.path.isdir(savedir):
    os.makedirs(savedir)


# Load protocols
protocol_0 = np.loadtxt(info.protocol_file_0)
protocol_0 = protocol_0.flatten().round()
protocol_0 = np.append([-80., 100.], protocol_0)  # Add 0.1s holding potential

protocol_1 = np.loadtxt(info.protocol_file_1)
protocol_1 = protocol_1.flatten().round()


# Create true model and synthetic data
model_true_0 = method.model.VCModel(info.true_model_file, transform=None, dt=dt, n_steps=n_steps)
model_true_1 = method.model.CCModel(info.true_model_file, transform=None, dt=dt, n_steps=n_steps)

parameters_true = np.ones(model_true_0.n_parameters())

model_true_0.design(protocol_0)
model_true_1.design(protocol_1)

times_0 = model_true_0.times()
data_0 = model_true_0.simulate(parameters_true, times_0)
data_0 += np.random.normal(0, info.noise_sigma_0, size=data_0.shape)
times_1 = model_true_1.times()
data_1 = model_true_1.simulate(parameters_true, times_1)
data_1 += np.random.normal(0, info.noise_sigma_1, size=data_1.shape)

if '--debug' in sys.argv:
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, figsize=(6,4))
    axes[0].plot(times_0, data_0)
    axes[0].set_ylabel('Current (A/F)')
    axes[0].set_xlabel('Time (ms)')
    axes[1].plot(times_1, data_1)
    axes[1].set_ylabel('Voltage (mV)')
    axes[1].set_xlabel('Time (ms)')
    plt.show()
    sys.exit()
del(model_true_0)
del(model_true_1)


# Estimate noise
noise_sigma_0 = np.std(data_0[:200])
print('Estimated noise sigma 0: ', noise_sigma_0)
noise_sigma_1 = np.std(data_1[:200])
print('Estimated noise sigma 1: ', noise_sigma_1)

# Create fitting model
model_fit_0 = method.model.VCModel(info.fit_model_file, transform=None, dt=dt, n_steps=n_steps)
model_fit_0.design(protocol_0)
model_fit_1 = method.model.CCModel(info.fit_model_file, transform=None, dt=dt, n_steps=n_steps)
model_fit_1.design(protocol_1)
if '--debug2' in sys.argv:
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, figsize=(6,4))
    axes[0].plot(times_0, data_0)
    axes[0].plot(times_0, model_fit_0.simulate(parameters_true))
    axes[0].set_ylabel('Current (A/F)')
    axes[0].set_xlabel('Time (ms)')
    axes[1].plot(times_1, data_1)
    axes[1].plot(times_1, model_fit_1.simulate(parameters_true))
    axes[1].set_ylabel('Voltage (mV)')
    axes[1].set_xlabel('Time (ms)')
    plt.show()
    sys.exit()

# Create likelihood
problem_0 = pints.SingleOutputProblem(model_fit_0, times_0, data_0)
log_likelihood_0 = pints.GaussianKnownSigmaLogLikelihood(problem_0, noise_sigma_0)
problem_1 = pints.SingleOutputProblem(model_fit_1, times_1, data_1)
log_likelihood_1 = pints.GaussianKnownSigmaLogLikelihood(problem_1, noise_sigma_1)
log_likelihood = pints.SumOfIndependentLogPDFs([log_likelihood_0, log_likelihood_1])
log_prior = pints.ComposedLogPrior(
    *([pints.LogNormalLogPrior(0, 0.15)] * problem_0.n_parameters())
)
transformation = pints.LogTransformation(problem_0.n_parameters())
log_posterior = pints.LogPosterior(log_likelihood, log_prior)


# Check log_posterior is not throwing error and deterministic
for _ in range(3):
    assert(log_posterior(np.ones(log_posterior.n_parameters())) ==\
            log_posterior(np.ones(log_posterior.n_parameters())))

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

ref_parameters = parameters_true

# Plot
pints.plot.pairwise(chains_final[0], kde=False, ref_parameters=ref_parameters)
plt.savefig('%s-fig1.png' % saveas)
plt.close('all')

pints.plot.trace(chains_final, ref_parameters=ref_parameters)
plt.savefig('%s-fig2.png' % saveas)
plt.close('all')

pints.plot.series(chains_final[0], problem)
plt.savefig('%s-fig3.png' % saveas)
plt.close('all')

# Check convergence using rhat criterion
print('R-hat:')
r_hats = pints.rhat(chains_final)
print(r_hats)
np.savetxt('%s-rhat.txt' % saveas, r_hats)
