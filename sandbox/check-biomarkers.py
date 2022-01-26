#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
sys.path.append('..')
import method.utils as utils
from method.model import CCModel
from method.biomarkers import BiomarkerExtractor

ds = 100  # downsampling factor when fitting splines

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

model.design([1000])

times = model.times()
aps = model.simulate(parameters, times)

# Try biomarkers
extractor = BiomarkerExtractor()
extractor.set_data(times, aps)
extractor.set_downsampling(ds)
biomarkers = extractor.available_biomarkers()
print('Available biomarkers: ', biomarkers)
extracted = extractor.extract()
for b in biomarkers:
    print(b, len(extracted[b]), extracted[b])

# Set up plot
plt.figure(figsize=(7, 4))
plt.plot(times, aps)
spl = utils.fit_spline(times, aps, ds=ds)
plt.plot(times, spl(times))
plt.plot(times, spl(times, 1))
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')

# Utils
baseline = utils.find_baseline(aps)
plt.axhline(baseline)

peak = utils.find_peak(aps)
plt.axhline(peak)

peak_idx_full = utils.find_peak_indices(times, aps, ds=ds)
plt.plot(times[peak_idx_full], aps[peak_idx_full], 'o')

base_idx_full = utils.find_baseline_indices(times, aps, ds=ds)
plt.plot(times[base_idx_full], aps[base_idx_full], 'o')

top_idx = utils.find_take_off_potential_indices(times, aps, ds=ds)
plt.plot(times[top_idx], aps[top_idx], 'o')

plt.show()
