#
# Biomarker functions for this project.
#
import numpy as np
import method.utils as utils


class BiomarkerExtractor(object):
    """
    A class for extracting biomarkers of spontaneous beating cardiac action
    potentials (APs).

    - rmp: (resting membrane potential) the most negative membrane potentials.
    - vmax: (action potential peak) the most positive membrane potentials.
    - mur: (max upstroke rate) the maximum values of the first derivative.
    - apa: (action potential amplitude) the voltage difference between the Vmax
           and the subsequent MDP.
    - apd50: the interval between the TOP and 50% repolarization.
    - apd90: the interval between the TOP and 90% repolarization.
    - tri: (triangulation ratio) difference between APD90 and APD50.

    More for non-adult-ventricular/atrial cells:
    - top: (take-off potential) the membrane potential when the first
           derivative of voltage with respect to time (dV/dt) reached 10% of
           its maximum value.
    - cl: (cycle length) the interval between MDPs in successive APs.
    - mrr: (max repolarisation rate) the minimum values of the first
           derivative.
    - trir: (triangulation ratio) ratio between APD50 and APD90.
           Defined in https://doi.org/10.1097/01.fjc.0000246853.15926.d4.
    - dd: (diastolic duration) the interval between MDP/RMP and TOP.
    - eddr: (early diastolic depolarisation rate) the slope of a linear fit
            between 10% and 50% of the diastolic duration.

    - early_diastolic_duration: the 10% and 50% of the diastolic duration.
    - late_diastolic_duration: the duration between 1% and 10% dV/dt.
    """

    BIOMARKER_LONG = dict(
        rmp = 'Resting membrane potential',
        vmax = 'Maximum potential peak',
        mur = 'Maximum upstroke rate',
        apa = 'Action potential amplitude',
        apd50 = r'APD$_{50}$',
        apd90 = r'APD$_{90}$',
        tri = 'Triangulation ratio',
    )

    BIOMARKER_SHORT = dict(
        rmp = 'RMP',
        vmax = r'$V_{max}$',
        mur = 'MUR',
        apa = 'APA',
        apd50 = r'APD$_{50}$',
        apd90 = r'APD$_{90}$',
        tri = 'Tri',
    )

    BIOMARKER_UNIT = dict(
        rmp = 'mV',
        vmax = 'mV',
        mur = 'V/s',
        apa = 'mV',
        apd50 = 'ms',
        apd90 = 'ms',
        tri = 'ms',
    )

    def __init__(self, pacing=None, times=None, aps=None, ds=None):
        self.BIOMARKERS = dict(
            rmp = self.rmp,
            vmax = self.vmax,
            mur = self.mur,
            apa = self.apa,
            apd50 = self.apd50,
            apd90 = self.apd90,
            tri = self.tri,
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

    def apd50(self, kwargs={}):
        """
        Extract the action potential durations at 50% repolarisation (APD50s).
        """
        return self.apdx(0.5, kwargs)

    def apd90(self, kwargs={}):
        """
        Extract the action potential durations at 90% repolarisation (APD50s).
        """
        return self.apdx(0.9, kwargs)

    def apdx(self, x, kwargs={}):
        kwargs = self._parse_kwargs_spline(kwargs)
        vx = self.aps[self._pidx] - x * self.apa(kwargs)
        apdx = []
        for i, (pi, fi) in enumerate(zip(self._pidx, self._fidx)):
            # idx = np.argmin(np.abs(self._spl(self.times[pi:fi]) - vx[i]))
            idx = np.argmin(np.abs(self.aps[pi:fi] - vx[i]))
            idx += pi
            if np.abs(self.aps[idx] - vx[i]) > 5:
                # Probably not hitting APDx (maybe altnan)
                apdx.append(np.nan)
            else:
                apdx.append(self.times[idx] - self.times[self._iidx[i]])
        return apdx

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

    def rmp(self, kwargs={}):
        """
        Extract resting membrane potential.
        """
        assert(self._iidx[0] > 0)
        return np.mean(self.aps[:self._iidx[0]])

    def tri(self, p1=0.9, p2=0.4, kwargs={}):
        """
        Extract triangulation of the action potentials.
        """
        tri = np.array(self.apdx(p1, kwargs)) - np.array(self.apdx(p2, kwargs))
        return tri

    def vmax(self, kwargs={}):
        """
        Extract peak values of the action potentials (V_{max}).
        """
        vmax = self.aps[self._pidx]
        return vmax
