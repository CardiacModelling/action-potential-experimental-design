#
# Utility functions
#
import pints
import numpy as np
from scipy.interpolate import UnivariateSpline

MIN_T = 180  # expected minimum duration between two APs.


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


def fit_spline(times, aps, ds=1, k=4, kwargs={}):
    # Fit and return a scipy.interpolate.UnivariateSpline object.
    #
    # times: (ms) array of time series.
    # aps: array of action potentials.
    # ds: downsampling factor when fitting splines.
    # k: degree of the spline.
    # kwargs: keyword arguments for scipy.interpolate.UnivariateSpline.
    #
    return UnivariateSpline(times[::ds], aps[::ds], k=k, **kwargs)


def find_baseline(aps, percentile=5):
    # Return an estimate of the baseline value for the action potentials.
    #
    # aps: array of action potentials.
    # percentile: percentile of the train of action potentials which defines
    #             the baseline.
    #
    return np.percentile(aps, percentile)


def find_baseline_indices(times, aps, spl=None, percentile=2, ds=1,
                          min_t=MIN_T):
    # Return the indices of the baseline for the action potentials.
    #
    # times: (ms) array of time series.
    # aps: array of action potentials.
    # spl: (scipy.interpolate.UnivariateSpline object) fitted spline of the
    #      data, refit if not provided.
    # percentile: percentile of the train of action potentials which defines
    #             the baseline.
    # ds: downsampling factor when fitting splines.
    # min_t: (ms) Duration which two peaks are not expected to occur.
    #
    dt = times[1] - times[0]
    min_idx = int((min_t) / dt)
    baseline = find_baseline(aps, percentile)
    vrange = 5  # mV
    if spl is None:
        spl = fit_spline(times, aps, ds=ds)
    r = spl.derivative(1).roots()
    idx = np.searchsorted(times, r, side='left')
    base_range = np.where((aps < baseline + vrange)
                          & (aps > baseline - vrange))[0]
    base_idx = np.intersect1d(idx, base_range)
    d_base_idx = base_idx[1:] - base_idx[:-1]
    base_idx = base_idx[np.append(True, d_base_idx > min_idx)]
    return base_idx


def find_peak(aps, percentile=95):
    # Return an estimate of the peak value for the action potentials.
    #
    # aps: array of action potentials.
    # percentile: percentile of the train of action potentials which defines
    #             the peak.
    #
    return np.percentile(aps, percentile)


def find_peak_indices(times, aps, spl=None, percentile=95, ds=1, min_t=MIN_T):
    # Return the indices of the peak for the action potentials.
    #
    # times: (ms) array of time series.
    # aps: array of action potentials.
    # spl: (scipy.interpolate.UnivariateSpline object) fitted spline of the
    #      data, refit if not provided.
    # percentile: percentile of the train of action potentials which defines
    #             the peak.
    # ds: downsampling factor when fitting splines.
    # min_t: (ms) Duration which two peaks are not expected to occur.
    #
    dt = times[1] - times[0]
    min_idx = int((min_t) / dt)
    peak = find_peak(aps, percentile)
    vrange = 15  # mV
    if spl is None:
        spl = fit_spline(times, aps, ds=ds)
    r = spl.derivative(1).roots()
    idx = np.searchsorted(times, r, side='left')
    peak_range = np.where((aps < peak + vrange) & (aps > peak - vrange))[0]
    peak_idx = np.intersect1d(idx, peak_range)
    d_peak_idx = peak_idx[1:] - peak_idx[:-1]
    peak_idx = peak_idx[np.append(True, d_peak_idx > min_idx)]
    return peak_idx


def find_take_off_potential_indices(times, aps, spl=None, method='slope', ds=1,
                                    min_t=MIN_T):
    # Return the indices of the take-off-potential for the action potentials.
    #
    # times: (ms) array of time series.
    # aps: array of action potentials.
    # spl: (scipy.interpolate.UnivariateSpline object) fitted spline of the
    #      data, refit if not provided.
    # method: ['intercept'|'slope'] method to find take-off-potential.
    # ds: downsampling factor when fitting splines.
    # min_t: (ms) Duration which two peaks are not expected to occur.
    #

    if spl is None:
        spl = fit_spline(times, aps, ds=ds)

    base_idx_full = find_baseline_indices(times, aps, spl=spl, ds=ds,
                                          min_t=min_t)
    peak_idx_full = find_peak_indices(times, aps, spl=spl, ds=ds, min_t=min_t)

    peak_idx = peak_idx_full.copy()
    base_idx = base_idx_full.copy()
    if peak_idx[0] < base_idx[0]:
        peak_idx = peak_idx_full[1:]
    if base_idx[-1] > peak_idx[-1]:
        base_idx = base_idx_full[:-1]
    #print(len(base_idx_full), len(peak_idx_full))
    #print(len(base_idx), len(peak_idx))
    assert(len(base_idx) == len(peak_idx))

    top_idx = []
    if method == 'intercept':
        for bi, pf in zip(base_idx, peak_idx):
            bf = bi + int(0.1 * (pf - bi))
            pi = bi + int(0.9 * (pf - bi))
            m1, b1 = np.polyfit(times[bi:bf], aps[bi:bf], 1)
            m2, b2 = np.polyfit(times[pi:pf], aps[pi:pf], 1)
            x = (b2 - b1) / (m1 - m2)
            top_idx.append(np.searchsorted(times, x, side='left'))
    elif method == 'slope':
        thres = 0.1  # following https://doi.org/10.1073/pnas.1308477110
        for bi, pf in zip(base_idx, peak_idx):
            t = times[bi:pf]
            dvdt = spl(t, 1)
            amaxdvdt = np.argmax(dvdt)
            # Find the closest value of the thres on the 'left' of the peak
            trimdvdt = dvdt[:amaxdvdt]
            ttop = t[:amaxdvdt][trimdvdt < thres * trimdvdt[-1]][-1]
            top_idx.append(np.searchsorted(times, ttop, side='left'))
    return top_idx


def get_n_action_potentials(times, aps, n=3, start=0, ds=1):
    # Return `n` full action potentials as tuples (times, aps).
    #
    # times: (ms) array of time series.
    # aps: array of action potentials.
    # n: number of full action potentials to return.
    # start: starting action potential.
    # ds: downsampling factor when fitting splines.
    #
    base_idx = find_baseline_indices(times, aps, ds=ds)
    start_idx = base_idx[start]
    end_idx = base_idx[start + n]
    return times[start_idx:end_idx], aps[start_idx:end_idx]


def find_first_mdp_index(times, aps, ds=1):
    # Return the first maximum diastolic potential (MDP) index.
    #
    # times: (ms) array of time series.
    # aps: array of action potentials.
    # ds: downsampling factor when fitting splines.
    #
    base_idx = find_baseline_indices(times, aps, ds=ds)
    return base_idx[0]


def get_n_action_potentials_by_peak(times, aps, n=3, start=0, ds=1):
    # Return `n` full action potentials as tuples (times, aps).
    #
    # times: (ms) array of time series.
    # aps: array of action potentials.
    # n: number of full action potentials to return.
    # start: starting action potential.
    # ds: downsampling factor when fitting splines.
    #
    peak_idx = find_peak_indices(times, aps, ds=ds)
    base_idx = find_baseline_indices(times, aps, ds=ds)
    if peak_idx[0] < base_idx[0]:
        peak_idx = peak_idx[1:]
    if base_idx[-1] > peak_idx[-1]:  # Though not necessary, still good check
        base_idx = base_idx[:-1]
    assert(len(base_idx) == len(peak_idx))
    start_idx = peak_idx[start]
    end_idx = base_idx[start + n]
    # end_idx = peak_idx[start + n]
    return times[start_idx:end_idx], aps[start_idx:end_idx]


def find_first_vmax_index(times, aps, ds=1):
    # Return the first action potential peak index.
    #
    # times: (ms) array of time series.
    # aps: array of action potentials.
    # ds: downsampling factor when fitting splines.
    #
    peak_idx = find_peak_indices(times, aps, ds=ds)
    return peak_idx[0]
