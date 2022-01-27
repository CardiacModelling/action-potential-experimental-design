#
# Utility functions
#
import pints
import numpy as np
from scipy.interpolate import UnivariateSpline


def save_protocol(fn, p, mode='vc'):
    # fn: file name to save as.
    # p: parameters to save.
    if mode == 'vc':
        with open(fn, 'w') as f:
            f.write('# Voltage [mV]\tDuration [ms]\n')
            for i in range(len(p) // 2):
                f.write(pints.strfloat(p[2 * i]) \
                        + '\t' \
                        + pints.strfloat(p[2 * i + 1]) \
                        + '\n' \
                        )
    elif mode == 'cc':
        with open(fn, 'w') as f:
            f.write('# Holding duration [ms]\n')
            for i in range(len(p)):
                f.write(pints.strfloat(p[i]) \
                        + '\n' \
                        )
    else:
        raise ValueError('mode must be either `vc` or `cc`.')


class PiecewiseSpline(object):
    def __init__(self, splines, intervals):
        assert(len(splines) == len(intervals))
        self._splines = splines
        self._intervals = np.asarray(intervals)

    def __call__(self, x, *args, **kwargs):
        # Look up for intervals
        #idx = ((self._intervals < x)[:, 0] & (self._intervals > x)[:, 1])
        #idx = np.where(idx)
        #return self._splines[idx](x, *args, **kwargs)
        out = np.full(len(x), np.nan)
        for i, s in zip(self._intervals, self._splines):
            idx = np.where((i[0] < x) & (i[1] >= x))[0]
            out[idx] = s(x[idx], *args, **kwargs)
        return out


def fit_spline(times, aps, pacing, ds=1, k=4, kwargs={}):
    # Fit piece-wise scipy.interpolate.UnivariateSpline to the action
    # potentials and return a PiecewiseSpline oject with the pacing indices,
    # action potential peak indices, and the next pacing indices (or the end).
    #
    # times: (ms) array of time series.
    # aps: array of action potentials.
    # ds: downsampling factor when fitting splines.
    # k: degree of the spline.
    # kwargs: keyword arguments for scipy.interpolate.UnivariateSpline.
    #
    spl = []
    interval = []
    iidx = []
    pidx = []
    fidx = []
    if pacing[0] > 0:
        idx_i = 0
        idx_f = np.argmin(np.abs(times - pacing[0]))
        interval.append((times[idx_i], times[idx_f]))
        spl.append(
            UnivariateSpline(
                times[idx_i:idx_f:ds],
                aps[idx_i:idx_f:ds],
                k=k,
                **kwargs
            )
        )
    for i in range(len(pacing)):
        idx_i = np.argmin(np.abs(times - pacing[i]))
        if i == len(pacing) - 1:
            idx_f = -1
        else:
            idx_f = np.argmin(np.abs(times - pacing[i + 1]))
        idx_p = idx_i + np.argmax(aps[idx_i:idx_f])
        interval.append((times[idx_i], times[idx_p]))
        spl.append(
            UnivariateSpline(
                times[idx_i:idx_p:ds],
                aps[idx_i:idx_p:ds],
                k=k,
                **kwargs
            )
        )
        interval.append((times[idx_p], times[idx_f]))
        spl.append(
            UnivariateSpline(
                times[idx_p:idx_f:ds],
                aps[idx_p:idx_f:ds],
                k=k,
                **kwargs
            )
        )
        iidx.append(idx_i)
        pidx.append(idx_p)
        fidx.append(idx_f)
    return PiecewiseSpline(spl, interval), iidx, pidx, fidx


def find_pacing_indices(times, aps, pacing):
    # Return the indices of pacing for the action potentials.
    #
    # times: (ms) array of time series.
    # aps: array of action potentials.
    # pacing: time at which pacing applied.
    #
    top_idx = []
    for p in pacing:
        top_idx.append(np.argmin(np.abs(times - p)))
    return top_idx
