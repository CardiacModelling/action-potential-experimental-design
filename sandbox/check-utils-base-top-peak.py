#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
sys.path.append('..')
import method.utils as utils
from method.model import CCModel

ds = 100  # downsampling factor when fitting splines

# Args
no_cal_top = '--no-top' in sys.argv
debug = '--debug' in sys.argv

#f = '../mmt/grandi-2010.mmt'
#f = '../mmt/ohara-2011.mmt'
#f = '../mmt/cipa-2017.mmt'
f = '../mmt/tomek-2019.mmt'
n_steps = 1
dt = 0.1

model = CCModel(
    f,
    transform=None,
    dt=dt,
    n_steps=n_steps,
)

parameters = [1] * model.n_parameters()

pacing = [50, 1000]
model.design(pacing[1:])

times = model.times()
aps = model.simulate(parameters, times)

# Set up plot
plt.figure(figsize=(7, 4))
plt.plot(times, aps)
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')

# Try utils
baseline = utils.find_baseline(aps)
plt.axhline(baseline)

peak = utils.find_peak(aps)
plt.axhline(peak)

if debug:
    import numpy as np
    from scipy.interpolate import UnivariateSpline
    spl = UnivariateSpline(times[::ds], aps[::ds], k=4)
    plt.plot(times, spl(times))
    plt.plot(times, spl(times, 1))
    r = spl.derivative(1).roots()
    idx = np.searchsorted(times, r, side='left')
    # plt.plot(r, np.zeros(len(r)), '.')

    peak_range = np.where((aps < peak + 5) & (aps > peak - 5))[0]
    peak_idx = np.intersect1d(idx, peak_range)
    plt.plot(times[peak_idx], aps[peak_idx], '.')

    base_range = np.where((aps < baseline + 5) & (aps > baseline - 5))[0]
    base_idx = np.intersect1d(idx, base_range)
    plt.plot(times[base_idx], aps[base_idx], '.')

peak_idx_full = utils.find_peak_indices(times, aps, ds=ds)
plt.plot(times[peak_idx_full], aps[peak_idx_full], 'o')

base_idx_full = utils.find_baseline_indices(times, aps, ds=ds)
plt.plot(times[base_idx_full], aps[base_idx_full], 'o')

if not no_cal_top:
    top_idx = utils.find_take_off_potential_indices(times, aps, pacing=pacing)
    print(top_idx)
    plt.plot(times[top_idx], aps[top_idx], 'o')

plt.show()
