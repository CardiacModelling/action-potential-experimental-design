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
    - apd: (action potential duration) the interval between the TOP and the
           subsequent RMP.
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
        mdp = 'Maximum diastolic potential',
        top = 'Take-off potential',
        vmax = 'Maximum potential peak',
        cl = 'Cycle length',
        mur = 'Maximum upstroke rate',
        mrr = 'Maximum repolarisation rate',
        apa = 'Action potential amplitude',
        apd = 'Action potential duration',
        apd50 = r'APD$_{50}$',
        apd90 = r'APD$_{90}$',
        tri = 'Triangulation ratio',
        dd = 'Diastolic duration',
        eddr = 'Early diastolic depolarisation rate',
    )

    BIOMARKER_SHORT = dict(
        mdp = 'MDP',
        top = 'TOP',
        vmax = r'$V_{max}$',
        cl = 'CL',
        mur = 'MUR',
        mrr = 'MRR',
        apa = 'APA',
        apd = 'APD',
        apd50 = r'APD$_{50}$',
        apd90 = r'APD$_{90}$',
        tri = 'Tri',
        dd = 'DD',
        eddr = 'EDDR',
    )

    BIOMARKER_UNIT = dict(
        mdp = 'mV',
        top = 'mV',
        vmax = 'mV',
        cl = 'ms',
        mur = 'V/s',
        mrr = 'V/s',
        apa = 'mV',
        apd = 'ms',
        apd50 = 'ms',
        apd90 = 'ms',
        tri = '1',
        dd = 'ms',
        eddr = 'V/s',
    )

    def __init__(self, times=None, aps=None, ds=None):
        self.BIOMARKERS = dict(
            mdp = self.mdp,
            top = self.top,
            vmax = self.vmax,
            cl = self.cl,
            mur = self.mur,
            mrr = self.mrr,
            apa = self.apa,
            apd = self.apd,
            apd50 = self.apd50,
            apd90 = self.apd90,
            tri = self.tri,
            dd = self.dd,
            eddr = self.eddr,
        )
        self._has_data = False
        self.set_downsampling(ds)
        self.set_one_biomarker_per_full_action_potential(True)
        self.set_data(times, aps)
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
            self._spl = utils.fit_spline(self.times, self.aps, ds=ds)
        else:
            self._spl = None

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

    def set_data(self, times, aps):
        """
        Set data for analysis.
        """
        self.times = times
        self.aps = aps
        self.has_data()
        self.fit_spline()
        if self._same_n_biomarkers and self._has_data:
            # detect number of full action potentials
            self._detect_n_aps()

    def set_downsampling(self, ds=None):
        """
        Set downsampling rate for fitting splines.
        """
        if ds is not None:
            self.ds = int(ds)
        else:
            self.ds = ds
        self.fit_spline()

    def set_one_biomarker_per_full_action_potential(self, option):
        """
        If True, extract biomarker only for a full action potential is
        detected.
        """
        self._same_n_biomarkers = bool(option)
        if self._same_n_biomarkers and self._has_data:
            # detect number of full action potentials
            self._detect_n_aps()

    def _apply_mask(self, x, m):
        """
        Apply filtering mask, which return x[m[0]:m[1]].
        """
        return x[m[0]:m[1]]

    def _detect_n_aps(self):
        """
        Detect number of full action potentials and masks for pidx and tidx.
        """
        kwargs = self._parse_kwargs_spline({})
        bidx = utils.find_baseline_indices(self.times, self.aps, **kwargs)
        tidx = utils.find_take_off_potential_indices(self.times, self.aps,
                                                     **kwargs)
        pidx = utils.find_peak_indices(self.times, self.aps, **kwargs)
        if tidx[0] < bidx[0]:
            ti = 1
        else:
            ti = 0
        if pidx[0] < bidx[0]:
            pi = 1
        else:
            pi = 0
        if tidx[-1] > bidx[-1]:
            tf = -1
        else:
            tf = None
        if pidx[-1] > bidx[-1]:
            pf = -1
        else:
            pf = None
        self._mask_tidx = (ti, tf)
        self._mask_pidx = (pi, pf)
        self._n_aps = len(bidx) - 1

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
        kwargs = self._parse_kwargs_spline(kwargs)
        pidx = utils.find_peak_indices(self.times, self.aps, **kwargs)
        bidx = utils.find_baseline_indices(self.times, self.aps, **kwargs)
        if self._same_n_biomarkers:
            pidx = self._apply_mask(pidx, self._mask_pidx)
        if bidx[0] < pidx[0]:
            bidx = bidx[1:]
        if pidx[-1] > bidx[-1]:
            pidx = pidx[:-1]
        assert(len(bidx) == len(pidx))
        apa = self.aps[pidx] - self.aps[bidx]
        if self._same_n_biomarkers:
            assert(self._n_aps == len(apa))
        return apa

    def apd(self, kwargs={}):
        """
        Extract the action potential durations (APDs).
        """
        kwargs = self._parse_kwargs_spline(kwargs)
        tidx = utils.find_take_off_potential_indices(self.times, self.aps,
                                                     **kwargs)
        bidx = utils.find_baseline_indices(self.times, self.aps, **kwargs)
        if self._same_n_biomarkers:
            tidx = self._apply_mask(tidx, self._mask_tidx)
        if bidx[0] < tidx[0]:
            bidx = bidx[1:]
        if tidx[-1] > bidx[-1]:
            tidx = tidx[:-1]
        assert(len(bidx) == len(tidx))
        apd = self.times[bidx] - self.times[tidx]
        if self._same_n_biomarkers:
            assert(self._n_aps == len(apd))
        return apd

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
        tidx = utils.find_take_off_potential_indices(self.times, self.aps,
                                                     **kwargs)
        bidx = utils.find_baseline_indices(self.times, self.aps, **kwargs)
        pidx = utils.find_peak_indices(self.times, self.aps, **kwargs)
        if self._same_n_biomarkers:
            tidx = self._apply_mask(tidx, self._mask_tidx)
            pidx = self._apply_mask(pidx, self._mask_pidx)
        if pidx[0] < tidx[0]:
            pidx = pidx[1:]
        if pidx[-1] > bidx[-1]:
            pidx = pidx[:-1]
        if bidx[0] < tidx[0]:
            bidx = bidx[1:]
        if tidx[-1] > bidx[-1]:
            tidx = tidx[:-1]
        assert(len(bidx) == len(tidx))
        assert(len(bidx) == len(pidx))
        vx = self.aps[pidx] - x * (self.aps[pidx] - self.aps[bidx])
        apdx = []
        for i, (pi, bi) in enumerate(zip(pidx, bidx)):
            idx = np.argmin(np.abs(self._spl(self.times[pi:bi]) - vx[i]))
            idx += pi
            apdx.append(self.times[idx] - self.times[tidx[i]])
        if self._same_n_biomarkers:
            assert(self._n_aps == len(apdx))
        return apdx

    def cl(self, kwargs={}):
        """
        Extract cycle lengths (CLs) of the action potentials.
        """
        kwargs = self._parse_kwargs_spline(kwargs)
        bidx = utils.find_baseline_indices(self.times, self.aps, **kwargs)
        mdp_times = self.times[bidx]
        cl = mdp_times[1:] - mdp_times[:-1]
        if self._same_n_biomarkers:
            assert(self._n_aps == len(cl))
        return cl

    def dd(self, kwargs={}):
        """
        Extract diastolic duration of the action potentials.
        """
        kwargs = self._parse_kwargs_spline(kwargs)
        bidx = utils.find_baseline_indices(self.times, self.aps, **kwargs)
        tidx = utils.find_take_off_potential_indices(self.times, self.aps,
                                                     **kwargs)
        if self._same_n_biomarkers:
            tidx = self._apply_mask(tidx, self._mask_tidx)
        if tidx[0] < bidx[0]:
            tidx = tidx[1:]
        if bidx[-1] > tidx[-1]:
            bidx = bidx[:-1]
        assert(len(bidx) == len(tidx))
        dd = self.times[tidx] - self.times[bidx]
        if self._same_n_biomarkers:
            assert(self._n_aps == len(dd))
        return dd

    def eddr(self, kwargs={}):
        """
        Extract early diastolic depolarisation rates of the action potentials.
        """
        kwargs = self._parse_kwargs_spline(kwargs)
        bidx = utils.find_baseline_indices(self.times, self.aps, **kwargs)
        tidx = utils.find_take_off_potential_indices(self.times, self.aps,
                                                     **kwargs)
        if self._same_n_biomarkers:
            tidx = self._apply_mask(tidx, self._mask_tidx)
        if tidx[0] < bidx[0]:
            tidx = tidx[1:]
        if bidx[-1] > tidx[-1]:
            bidx = bidx[:-1]
        assert(len(bidx) == len(tidx))
        eddr = []
        for b, t in zip(bidx, tidx):
            # Get 10% to 50% of diastolic durations
            times = self.times[b:t]
            voltage = self.aps[b:t]
            i10 = int(0.1 * len(times))
            i50 = int(0.5 * len(times))
            # Fit a line
            m, b = np.polyfit(times[i10:i50], voltage[i10:i50], 1)
            eddr.append(m)
        if self._same_n_biomarkers:
            assert(self._n_aps == len(eddr))
        return eddr

    def mrr(self, kwargs={}):
        """
        Extract maximum rates of the action potential repolarisations.
        """
        kwargs = self._parse_kwargs_spline(kwargs)
        bidx = utils.find_baseline_indices(self.times, self.aps, **kwargs)
        pidx = utils.find_peak_indices(self.times, self.aps, **kwargs)
        if self._same_n_biomarkers:
            pidx = self._apply_mask(pidx, self._mask_pidx)
        if bidx[0] < pidx[0]:
            bidx = bidx[1:]
        if pidx[-1] > bidx[-1]:
            pidx = pidx[:-1]
        assert(len(bidx) == len(pidx))
        mrr = []
        for i, (b, p) in enumerate(zip(bidx, pidx)):
            t = self.times[p:b]
            mrr.append(np.min(self._spl(t, 1)))  # 1 means first derivative
        if self._same_n_biomarkers:
            assert(self._n_aps == len(mrr))
        return mrr

    def mur(self, kwargs={}):
        """
        Extract maximum rates of the action potential upstrokes.
        """
        kwargs = self._parse_kwargs_spline(kwargs)
        bidx = utils.find_baseline_indices(self.times, self.aps, **kwargs)
        pidx = utils.find_peak_indices(self.times, self.aps, **kwargs)
        if self._same_n_biomarkers:
            pidx = self._apply_mask(pidx, self._mask_pidx)
        if pidx[0] < bidx[0]:
            pidx = pidx[1:]
        if bidx[-1] > pidx[-1]:
            bidx = bidx[:-1]
        assert(len(bidx) == len(pidx))
        mur = []
        for i, (b, p) in enumerate(zip(bidx, pidx)):
            t = self.times[b:p]
            mur.append(np.max(self._spl(t, 1)))  # 1 means first derivative
        if self._same_n_biomarkers:
            assert(self._n_aps == len(mur))
        return mur

    def mdp(self, kwargs={}):
        """
        Extract maximum diastolic potentials (MDPs) of the action potentials.
        """
        kwargs = self._parse_kwargs_spline(kwargs)
        bidx = utils.find_baseline_indices(self.times, self.aps, **kwargs)
        mdp = self.aps[bidx[1:]]
        if self._same_n_biomarkers:
            assert(self._n_aps == len(mdp))
        return mdp

    def top(self, kwargs={}):
        """
        Extract take off potentials (TOPs) of the action potentials.
        """
        kwargs = self._parse_kwargs_spline(kwargs)
        tidx = utils.find_take_off_potential_indices(self.times, self.aps,
                                                     **kwargs)
        if self._same_n_biomarkers:
            tidx = self._apply_mask(tidx, self._mask_tidx)
        top = self.aps[tidx]
        if self._same_n_biomarkers:
            assert(self._n_aps == len(top))
        return top

    def tri(self, kwargs={}):
        """
        Extract triangulation ratio of the action potentials.
        """
        tri = np.array(self.apd50(kwargs)) / np.array(self.apd90(kwargs))
        if self._same_n_biomarkers:
            assert(self._n_aps == len(tri))
        return tri

    def vmax(self, kwargs={}):
        """
        Extract peak values of the action potentials (V_{max}).
        """
        kwargs = self._parse_kwargs_spline(kwargs)
        pidx = utils.find_peak_indices(self.times, self.aps, **kwargs)
        if self._same_n_biomarkers:
            pidx = self._apply_mask(pidx, self._mask_pidx)
        vmax = self.aps[pidx]
        if self._same_n_biomarkers:
            assert(self._n_aps == len(vmax))
        return vmax
