#
# Biomarker functions for this project.
#
import numpy as np
import method.utils as utils


class BiomarkerExtractor(object):
    """
    A class for extracting biomarkers of spontaneous beating cardiac action
    potentials (APs).

    - apa: (action potential amplitude) the voltage difference between pacing
           and Vmax.
    - apd10: the interval between pacing and 10% repolarization.
    - apd20: the interval between pacing and 20% repolarization.
    - apd30: the interval between pacing and 30% repolarization.
    - apd40: the interval between pacing and 40% repolarization.
    - apd50: the interval between pacing and 50% repolarization.
    - apd60: the interval between pacing and 60% repolarization.
    - apd70: the interval between pacing and 70% repolarization.
    - apd80: the interval between pacing and 80% repolarization.
    - apd90: the interval between pacing and 90% repolarization.
    - dp: (dome peak) the most positive value during pd.
    - mur: (max upstroke rate) the maximum values of the first derivative.
    - pd: (plateau duration) duration at which |first derivative| < 0.15 V/s.
    - rmp: (resting membrane potential) the most negative membrane potentials.
    - t2p: (time to peak) duration between pacing and vmax.
    - tri: (triangulation) difference between APD90 and APD40.
    - varea: (area under action potentials) Area under the action potential
             within the period defined by APD30.
    - vmax: (action potential peak) the most positive membrane potentials.

    More for non-adult-ventricular/atrial cells:
    - cl: (cycle length) the interval between MDPs in successive APs.
    - dd: (diastolic duration) the interval between MDP/RMP and TOP.
    - eddr: (early diastolic depolarisation rate) the slope of a linear fit
            between 10% and 50% of the diastolic duration.
    - early_diastolic_duration: the 10% and 50% of the diastolic duration.
    - late_diastolic_duration: the duration between 1% and 10% dV/dt.
    - mrr: (max repolarisation rate) the minimum values of the first
           derivative.
    - top: (take-off potential) the membrane potential when the first
           derivative of voltage with respect to time (dV/dt) reached 10% of
           its maximum value.
    - trir: (triangulation ratio) ratio between APD50 and APD90.
           Defined in https://doi.org/10.1097/01.fjc.0000246853.15926.d4.

    """

    BIOMARKER_LONG = dict(
        apa = 'Action potential amplitude',
        apd10 = r'APD$_{10}$',
        apd20 = r'APD$_{20}$',
        apd30 = r'APD$_{30}$',
        apd40 = r'APD$_{40}$',
        apd50 = r'APD$_{50}$',
        apd60 = r'APD$_{60}$',
        apd70 = r'APD$_{70}$',
        apd80 = r'APD$_{80}$',
        apd90 = r'APD$_{90}$',
        dp = 'Dome peak',
        mur = 'Maximum upstroke rate',
        pd = 'Plateau duration',
        rmp = 'Resting membrane potential',
        t2p = 'Time to peak',
        tri = 'Triangulation ratio',
        varea = 'Area under action potential',
        vmax = 'Maximum potential peak',
    )

    BIOMARKER_SHORT = dict(
        apa = 'APA',
        apd10 = r'APD$_{10}$',
        apd20 = r'APD$_{20}$',
        apd30 = r'APD$_{30}$',
        apd40 = r'APD$_{40}$',
        apd50 = r'APD$_{50}$',
        apd60 = r'APD$_{60}$',
        apd70 = r'APD$_{70}$',
        apd80 = r'APD$_{80}$',
        apd90 = r'APD$_{90}$',
        dp = 'DP',
        mur = 'MUR',
        pd = 'PD',
        rmp = 'RMP',
        t2p = 'T2P',
        tri = 'Tri',
        varea = r'$V_{area}$',
        vmax = r'$V_{max}$',
    )

    BIOMARKER_UNIT = dict(
        apa = 'mV',
        apd10 = 'ms',
        apd20 = 'ms',
        apd30 = 'ms',
        apd40 = 'ms',
        apd50 = 'ms',
        apd60 = 'ms',
        apd70 = 'ms',
        apd80 = 'ms',
        apd90 = 'ms',
        dp = 'mV',
        mur = 'V/s',
        pd = 'ms',
        rmp = 'mV',
        t2p = 'ms',
        tri = 'ms',
        varea = 'mVms',
        vmax = 'mV',
    )

    def __init__(self, pacing=None, times=None, aps=None, ds=None):
        self.BIOMARKERS = dict(
            apa = self.apa,
            apd10 = self.apd10,
            apd20 = self.apd20,
            apd30 = self.apd30,
            apd40 = self.apd40,
            apd50 = self.apd50,
            apd60 = self.apd60,
            apd70 = self.apd70,
            apd80 = self.apd80,
            apd90 = self.apd90,
            dp = self.dp,
            mur = self.mur,
            pd = self.pd,
            rmp = self.rmp,
            t2p = self.t2p,
            tri = self.tri,
            varea = self.varea,
            vmax = self.vmax,
        )
        self._has_data = False
        self.set_downsampling(ds)
        self.set_data(pacing, times, aps)
        self.set_biomarker_list(self.available_biomarkers())

    def available_biomarkers(self):
        """
        Return a list of available biomarkers.
        """
        return list(self.BIOMARKERS.keys())

    def biomarker_long_name(self, b):
        """
        Return the long name of the biomarker b.
        """
        return self.BIOMARKER_LONG[b]

    def biomarker_short_name(self, b):
        """
        Return the short name of the biomarker b.
        """
        return self.BIOMARKER_SHORT[b]

    def biomarker_unit(self, b):
        """
        Return the unit of the biomarker b.
        """
        return self.BIOMARKER_UNIT[b]

    def extract(self):
        """
        Extract biomarkers of the action potentials.
        """
        self.has_data()
        if not self._has_data:
            raise ValueError('Incorrect data format or no data is provided.')

        d = {}
        for k in self.list_of_biomarkers:
            v = self.BIOMARKERS[k]
            d[k] = v()
        return d

    def fit_spline(self):
        if self._has_data:
            if self.ds is not None:
                ds = self.ds
            else:
                ds = 1
            self._spl, self._iidx, self._pidx, self._fidx = utils.fit_spline(
                self.times,
                self.aps,
                pacing=self.pacing,
                ds=ds,
            )
        else:
            self._spl = None
            self._iidx = None
            self._pidx = None
            self._fidx = None

    def has_data(self):
        """
        Simple check that data are acceptable input.
        """
        try:
            if (len(self.times) == len(self.aps)):
                self._has_data = True
            else:
                self._has_data = False
        except TypeError:
            self._has_data = False

    def set_biomarker_list(self, list_of_biomarkers):
        self.list_of_biomarkers = list_of_biomarkers

    def set_data(self, pacing, times, aps):
        """
        Set data for analysis.
        """
        self.pacing = pacing
        self.times = times
        self.aps = aps
        self.has_data()
        self.fit_spline()

    def set_downsampling(self, ds=None):
        """
        Set downsampling rate for fitting splines.
        """
        if ds is not None:
            self.ds = int(ds)
        else:
            self.ds = ds
        self.fit_spline()

    def _parse_kwargs_spline(self, kwargs):
        if self.ds is not None:
            kwargs['ds'] = self.ds
        if self._spl is not None:
            kwargs['spl'] = self._spl
        return kwargs

    def apa(self, kwargs={}):
        """
        Extract the action potential heights.
        """
        apa = self.aps[self._pidx] - self.aps[self._iidx]
        return apa

    def apd10(self, kwargs={}):
        """
        Extract the action potential durations at 10% repolarisation (APD50s).
        """
        return self.apdx(0.1, kwargs)

    def apd20(self, kwargs={}):
        """
        Extract the action potential durations at 20% repolarisation (APD50s).
        """
        return self.apdx(0.2, kwargs)

    def apd30(self, kwargs={}):
        """
        Extract the action potential durations at 30% repolarisation (APD50s).
        """
        return self.apdx(0.3, kwargs)

    def apd40(self, kwargs={}):
        """
        Extract the action potential durations at 40% repolarisation (APD50s).
        """
        return self.apdx(0.4, kwargs)

    def apd50(self, kwargs={}):
        """
        Extract the action potential durations at 50% repolarisation (APD50s).
        """
        return self.apdx(0.5, kwargs)

    def apd60(self, kwargs={}):
        """
        Extract the action potential durations at 60% repolarisation (APD50s).
        """
        return self.apdx(0.6, kwargs)

    def apd70(self, kwargs={}):
        """
        Extract the action potential durations at 70% repolarisation (APD50s).
        """
        return self.apdx(0.7, kwargs)

    def apd80(self, kwargs={}):
        """
        Extract the action potential durations at 80% repolarisation (APD50s).
        """
        return self.apdx(0.8, kwargs)

    def apd90(self, kwargs={}):
        """
        Extract the action potential durations at 90% repolarisation (APD50s).
        """
        return self.apdx(0.9, kwargs)

    def apdx(self, x, kwargs={}):
        """
        Extract the action potential durations at x (fraction) repolarisation.
        """
        vx = self.aps[self._pidx] - x * self.apa(kwargs)
        apdx = []
        for i, (ii, pi, fi) in enumerate(zip(self._iidx, self._pidx, self._fidx)):
            t_ip = np.linspace(self.times[ii], self.times[pi], 100)
            idx_i = np.argmin(np.abs(self._spl(t_ip) - vx[i]))
            # idx_i = np.argmin(np.abs(self.aps[ii:pi] - vx[i]))
            idx_i += ii
            # idx_f = np.argmin(np.abs(self._spl(self.times[pi:fi]) - vx[i]))
            idx_f = np.argmin(np.abs(self.aps[pi:fi] - vx[i]))
            idx_f += pi
            if np.abs(self.aps[idx_f] - vx[i]) > 5:
                # Probably not hitting APDx (maybe altnan)
                apdx.append(np.nan)
            else:
                apdx.append(self.times[idx_f] - self.times[idx_i])
        return apdx

    def dp(self, thr=0.15, tmax=300, kwargs={}):
        """
        Extract dome peak values.
        """
        dp = []
        for pi in self._pidx:
            fi = np.argmin(np.abs(self.times - (self.times[pi] + tmax)))
            t = self.times[pi:fi]
            idx = np.where(np.abs(self._spl(t, 1)) < thr)[0]  # 1 means first derivative
            dp.append(np.max(self.aps[idx]))
        return dp

    def mrr(self, kwargs={}):
        """
        Extract maximum rates of the action potential repolarisations.
        """
        mrr = []
        for i, (pi, fi) in enumerate(zip(self._pidx, self._fidx)):
            t = self.times[pi:fi]
            mur.append(np.min(self._spl(t, 1)))  # 1 means first derivative
        return mrr

    def mur(self, kwargs={}):
        """
        Extract maximum rates of the action potential upstrokes.
        """
        mur = []
        for i, (ii, pi) in enumerate(zip(self._iidx, self._pidx)):
            t = self.times[ii:pi]
            mur.append(np.max(self._spl(t, 1)))  # 1 means first derivative
        return mur

    def pd(self, thr=0.15, tmax=300, kwargs={}):
        """
        Extract plateau durations.
        """
        pd = []
        for pi in self._pidx:
            fi = np.argmin(np.abs(self.times - (self.times[pi] + tmax)))
            t = self.times[pi:fi]
            idx = np.where(np.abs(self._spl(t, 1)) < thr)[0]  # 1 means first derivative
            dt = self.times[idx[1:]] - self.times[idx[:-1]]
            pd.append(np.sum(dt[np.abs(dt - np.min(dt)) < 1e-4]))
        return pd

    def rmp(self, kwargs={}):
        """
        Extract resting membrane potential.
        """
        assert(self._iidx[0] > 0)
        return np.mean(self.aps[:self._iidx[0]])

    def t2p(self, kwargs={}):
        """
        Extract time to peak values.
        """
        t2p = self.times[self._pidx] - self.times[self._iidx]
        return t2p

    def tri(self, p1=0.9, p2=0.4, kwargs={}):
        """
        Extract triangulation of the action potentials.
        """
        tri = np.array(self.apdx(p1, kwargs)) - np.array(self.apdx(p2, kwargs))
        return tri

    def varea(self, x=0.3, kwargs={}):
        """
        Compute the area under the action potentials.
        """
        vx = self.aps[self._pidx] - x * self.apa(kwargs)
        varea = []
        for i, (ii, pi, fi) in enumerate(zip(self._iidx, self._pidx, self._fidx)):
            t_ip = np.linspace(self.times[ii], self.times[pi], 100)
            idx_i = np.argmin(np.abs(self._spl(t_ip) - vx[i]))
            idx_i += ii
            # idx_f = np.argmin(np.abs(self._spl(self.times[pi:fi]) - vx[i]))
            idx_f = np.argmin(np.abs(self.aps[pi:fi] - vx[i]))
            idx_f += pi
            if np.abs(self.aps[idx_f] - vx[i]) > 5:
                # Probably not hitting APDx (maybe altnan)
                varea.append(np.nan)
            else:
                varea.append(np.trapz(self.aps[idx_i:idx_f],
                                      self.times[idx_i:idx_f]))
        return varea

    def vmax(self, kwargs={}):
        """
        Extract peak values of the action potentials (V_{max}).
        """
        vmax = self.aps[self._pidx]
        return vmax
