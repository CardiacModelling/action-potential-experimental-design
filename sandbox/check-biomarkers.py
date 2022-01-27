#!/usr/bin/env python3
import sys
sys.path.append('..')
import matplotlib.pyplot as plt
import numpy as np
import method.utils as utils
from method.model import CCModel
from method.biomarkers import BiomarkerExtractor

ds = 1  # downsampling factor when fitting splines

#f = '../mmt/grandi-2010.mmt'
#f = '../mmt/ohara-2011.mmt'
#f = '../mmt/cipa-2017.mmt'
f = '../mmt/tomek-2019.mmt'
n_steps = 0
dt = 0.1

model = CCModel(
    f,
    transform=None,
    dt=dt,
    n_steps=n_steps,
)

parameters = [1] * model.n_parameters()

pacing = [50]
model.design(pacing[1:])
pacing = np.cumsum(pacing)

times = model.times()
aps = model.simulate(parameters, times)

# Try biomarkers
extractor = BiomarkerExtractor()
extractor.set_data(pacing, times, aps)
extractor.set_downsampling(ds)
biomarkers = extractor.available_biomarkers()
print('Available biomarkers: ', biomarkers)
extracted = extractor.extract()
for b in biomarkers:
    print(b, extracted[b])

# Set up plot
plt.figure(figsize=(7, 4))
plt.plot(times, aps)
spl, _, _, _ = utils.fit_spline(times, aps, pacing, ds=ds)
plt.plot(times, spl(times))
plt.plot(times, spl(times, 1))
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')

plt.show()
