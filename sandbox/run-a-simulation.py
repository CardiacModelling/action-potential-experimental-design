#!/usr/bin/env python3
import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
import method.model

f = '../mmt/grandi-2010.mmt'
#f = '../mmt/ohara-2011.mmt'
#f = '../mmt/tomek-2019.mmt'
n_steps = 20
dt = 0.1

model = method.model.VCModel(
    f,
    transform=None,
    dt=dt,
    n_steps=n_steps,
)

parameters = np.ones(model.n_parameters())

opt_measure = 'LSA-A'
opt_model = 'ohara-2011'
pf = '../design/out/' + opt_measure + '-vc-' + opt_model \
    + '/opt-prt-run0-rank0.txt'
p = np.loadtxt(pf).flatten().round()
times = np.arange(0, np.sum(p[1::2]), dt)

model.set_voltage_protocol(p)

plt.plot(times, model.voltage(times))
plt.plot(times, model.simulate(parameters, times))
plt.show()
