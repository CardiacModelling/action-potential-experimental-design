#
# OED score function with simple Shannon entropy using GSA approximation.
#
import numpy as np
import pyoed


def shannon_all_measure(s):
    # Compute a Shannon entropy measure over the whole time horizon using a
    # given sensitivity matrix `s` (expecting GSA).
    # Shannon_all = 1. / sum_j(-1 * mean_i(s) ln(mean_i(s)))
    #
    # s: the sensitivity matrix (s_yi_thetaj)
    #    [s_y1_theta1, s_y1_theta2, ...]
    #    [s_y2_theta1, s_y2_theta2, ...]
    #    [...        , ...        , ...].
    ss = np.copy(s)
    m = np.mean(ss, axis=0)
    # Calculate p * log(p)
    # NOTE: although Sobol indices in theory larger than zero, we are
    #       estimating it using a finite (and potentially rather small) sample
    #       size. Therefore it's possible to get zero or even negative values
    #       due to the uncertainty.
    l0 = (m <= 0) | (m >= 1)
    m[l0] = 0  # Assume 0 * log(0) = 0 and 1 * log(1) = 0
    m[~l0] = -1. * m[~l0] * np.log(m[~l0])  # p * log(p)
    return 1. / np.sum(m, axis=0)


def shannon_segments_measure(s, t, seg):
    # Compute a Shannon entropy measure over a some time segments using a
    # given sensitivity matrix `s` (expecting GSA).
    # Shannon_seg = mean_i sum_j(-1 * mean_t_seg(s) ln(mean_t_seg(s)))
    #
    # s: the sensitivity matrix (s_yi_thetaj)
    #    [s_y1_theta1, s_y1_theta2, ...]
    #    [s_y2_theta1, s_y2_theta2, ...]
    #    [...        , ...        , ...].
    # t: an array of time points ti corresponding to each s_yi.
    # seg: time for each segment duration, assuming sum(seg) = last entry of t.
    #return np.mean(np.sum(s * np.log(s), axis=1), axis=0)  # each time point
    ss = np.copy(s)
    t_f = np.cumsum(seg)
    t_i = np.roll(t_f, 1)
    t_i[0] = 0  # assume time start from 0
    seg_idx = (t_i.reshape((-1, 1)) <= t) & (t < t_f.reshape((-1, 1)))
    m_seg = [np.mean(ss[seg_i, :], axis=0) for seg_i in seg_idx]
    m_seg = np.asarray(m_seg, dtype=np.float)
    # Calculate p * log(p)
    # NOTE: although Sobol indices in theory larger than zero, we are
    #       estimating it using a finite (and potentially rather small) sample
    #       size. Therefore it's possible to get zero or even negative values
    #       due to the uncertainty.
    l0 = (m_seg <= 0) | (m_seg >= 1)
    m_seg[l0] = 0  # Assume 0 * log(0) = 0 and 1 * log(1) = 0
    m_seg[~l0] = -1. * m_seg[~l0] * np.log(m_seg[~l0])  # p * log(p)
    return np.mean(np.sum(m_seg, axis=1), axis=0)


def parameter_dependency_measure(s):
    # Compute the parmeter dependency using a given sensitivity matrix `s`.
    # parameter_dependency = sum_i(1 - sum_j(s))
    #
    # s: the sensitivity matrix (s_yi_thetaj)
    #    [s_y1_theta1, s_y1_theta2, ...]
    #    [s_y2_theta1, s_y2_theta2, ...]
    #    [...        , ...        , ...].
    ss = np.copy(s)
    return np.sum(1. - np.sum(ss, axis=1), axis=0)


def output_uncertainty_measure(y):
    # Compute the overal output uncertainty using a given set of model output
    # `y`.
    # output_uncertainty = 1. / sum_i(var(y_i))
    #
    # y: set of model outputs (y_ti_thetaj)
    #    [y_t1_theta1, y_t1_theta2, ...]
    #    [y_t2_theta1, y_t2_theta2, ...]
    #    [...        , ...        , ...]
    yy = np.copy(y)
    return np.abs(1. / np.sum(np.var(yy, ddof=1, axis=1), axis=0))


class ShannonDesignMeasure(pyoed.DesignMeasure):
    """
    Self define error measure for square wave protocol optimisation.

    This tries to minimise the Shannon entropy via estimated using GSA sobol
    indice [1].
    
    [1] Schenkendorf, R., Xie, X., Rehbein, M., Scholl, S., and Krewer, U.
        (2018). The impact of global sensitivities and design measures in
        model-based optimal experimental design. Processes, 6(4):27.
    """
    
    def __init__(self, model, boundaries, method, method_kw={},
                 weight=[1, 1, 0, 1]):
        """
        # model: a pyoed.ForwardModel object.
        # boundaries: Boundaries of the model parameter space.
        # method: A pyoed.Sensitivity class.
        # method_kw: A dictionary with keyword arguments passed to the method
        #            constructor.
        # weight: Shannon-entropy based multi-objective design weighting
        #         vector for
        #         [shannon_all_measure,
        #          shannon_segments_measure,
        #          parameter_dependency_measure,
        #          output_uncertainty_measure].
        """
        super(ShannonDesignMeasure, self).__init__()

        self._model = model
        self._boundaries = boundaries
        self._n_mp = model.n_model_parameters()
        self._n_dp = model.n_design_variables()
        self.set_maximum_error(float('inf'))
        self.set_n_batches(None)
        if len(weight) != 4:
            raise ValueError('Weight must be a list of length 4.')
        self._weight = np.array(weight)

        # Create sensitivity method
        if method is not None and not issubclass(method, pyoed.Sensitivity):
            raise ValueError('Method must be subclass of pyoed.Sensitivity.')
        if not method:
            method = pyoed.SobolFirstOrderSensitivity
        self._method = method(self._boundaries, **method_kw)

    def set_n_batches(self, n):
        """
        Set number of batches for calculating the sensitivity. ``n`` must be an
        integer larger than 1 or ``None``. Default ``None``.
        """
        self._n_batches = n

    def n_batches(self):
        """
        Return the currently set number of batches for calculating the
        sensitivity.
        """
        return self._n_batches

    def set_maximum_error(self, e):
        """
        Set the maximum error value during ``__call__``, when the result is not
        finite.
        """
        self._max_error = e

    def maximum_error(self):
        """
        Return the currently set maximum error value during ``__call__``.
        """
        return self._max_error
 
    def n_parameters(self):
        return self._n_dp

    def sub_measures(self, param):
        # Return each sub Shannon measure

        # Get step duration
        seg_t = param[1::2]

        # Update design
        self._model.design(param)
        # Check if the protocol gives nonsense simulations
        # if not np.all(np.isfinite(self._model.simulated_currents)):
        #     return float('inf')

        # ask()
        ps = self._method.ask()
        fs = []
        for p in ps:
            fs.append(self._model.simulate(p))
        fs = np.asarray(fs)
        if not np.isfinite(fs).all():
            return self._max_error

        # tell()
        if self._n_batches is None:
            s = self._method.tell(fs)
        else:
            n_outputs = len(fs[0])
            s = np.zeros((self._n_mp, n_outputs))
            s[:] = np.nan
            n_each = n_outputs // self._n_batches
            start = 0
            if n_each > 0:
                for i in range(self._n_batches - 1):
                    end = (i + 1) * n_each
                    s[:, start:end] = self._method.tell(fs[:, start:end])
                    start = end
            s[:, start:] = self._method.tell(fs[:, start:])
            s = s.T  # rows being model outputs; cols being model parameters

        # Measure
        ds_t = self._model.times()  # time point for each s

        s1 = shannon_all_measure(s)
        s2 = shannon_segments_measure(s, ds_t, seg_t)
        s3 = parameter_dependency_measure(s)
        s4 = output_uncertainty_measure(fs.T)
        s = np.array([s1, s2, s3, s4]) * self._weight
        return s

    def __call__(self, param):
        # Return the design measure value.

        try:
            s = self.sub_measures(param)
            s = np.sum(s)
            if np.isfinite(s):
                return s
            else:
                return self._max_error
        except:
            return self._max_error  # Just make sure it won't break the opt.
