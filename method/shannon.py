#
# OED score function with simple Shannon entropy using GSA approximation.
#
from __future__ import print_function
import numpy as np
import pyoed
from SALib.sample import saltelli
import sobol


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
        if method is not None and not issubclass(method, pyoed.Sensitivity):
            raise ValueError('Method must be subclass of pyoed.Sensitivity.')
        self._method = method
        self.set_n_batches(None)
        self._n_samples = n_samples
        if len(weight) != 4:
            raise ValueError('Weight must be a list of length 4.')
        self._weight = np.array(weight)
 
    def n_parameters(self):
        return self._n_working_parameters

    def set_voltage(self, v=None):
        self._set_voltage = v
        if v is None and self._set_duration is None:
            self._n_working_parameters = self._n_parameters
        elif v is not None and self._set_duration is None:
            self._n_working_parameters = int(self._n_parameters / 2)
        else:
            raise ValueError('Both voltage and duration are fixed.')

    def set_duration(self, v=None):
        self._set_duration = v
        if v is None and self._set_voltage is None:
            self._n_working_parameters = self._n_parameters
        elif v is not None and self._set_voltage is None:
            self._n_working_parameters = int(self._n_parameters / 2)
        else:
            raise ValueError('Both voltage and duration are fixed.')

    def parameter_samples(self):
        # Return the parameter samples for calculating the Sobol sensitivity
        return self._parameter_samples

    def _convert_parameter(self, p):
        if self._set_voltage is not None and self._set_duration is not None:
            raise ValueError('Both voltage and duration are fixed.')

        if self._set_voltage is not None:
            param_duration = np.copy(p)
            p = np.zeros(len(param_duration))
            p[::2] = np.copy(self._set_voltage)
            p[1::2] = param_duration

        if self._set_duration is not None:
            param_voltage = np.copy(p)
            p = np.zeros(len(param_duration))
            p[1::2] = np.copy(self._set_duration)
            p[::2] = param_voltage

        return p

    def sub_measures(self, param, debug=False):
        # Return each sub Shannon measure

        # Get parameter
        param = self._convert_parameter(param)
        seg_t = param[1::2]  # step duration

        # Update protocol
        for model in self._model_list:
            model.set_voltage_protocol(param)
            # Check if the protocol gives nonsense simulations
            if not np.all(np.isfinite(model.simulated_currents)):
                return np.inf

        ds = 100  # NOTE downsample rate for the output time points

        scores = []
        for model in self._model_list:
            # Run simulations
            ## Run simulation per parameter sample
            #sims = []
            #for i, p in enumerate(self._parameter_samples):
            #    sims.append(model.simulate(p))
            #sims = np.asarray(sims)
            # Broadcast method to run simulations for all parameter samples
            sims = model.simulate(self._parameter_samples, downsample=ds,
                    multi_input=True)

            if not np.all(np.isfinite(sims)):
                return np.inf

            ds_time_points = len(model.times()[::ds])
            S_gsa = np.zeros((ds_time_points, self._n_model_parameters))
            ## Use SALib.analyze.sobol
            #for i, idx in enumerate(time_idx):
            #    # Compute Sobol sensitivity
            #    Si = sobol.analyze(self._salib_problem, sims[:, idx],
            #            calc_second_order=False)
            #    # Schenkendorf et al. 2018 (Eq. 10) uses only S1
            #    S_gsa[i, :] = Si['S1']
            # Use reimplemented sobol
            Si = sobol.analyze(self._salib_problem, sims, conf_level=None,
                    calc_second_order=False)

            S_gsa[:, :] = Si['S1'].T

            ds_t = model.times()[::100]  # time point for each S_gsa

            s1 = shannon_all_measure(S_gsa)
            s2 = shannon_segments_measure(S_gsa, ds_t, seg_t)
            s3 = parameter_dependency_measure(S_gsa)
            s4 = output_uncertainty_measure(sims.T)
            s = np.array([s1, s2, s3, s4]) * self._weight

            scores.append(s)

        # Take mean over all input models
        s = np.mean(scores, axis=0)
        return s

    def __call__(self, param, debug=False):

        try:
            # Get parameter
            param = self._convert_parameter(param)
            seg_t = param[1::2]  # step duration

            # Update protocol
            for model in self._model_list:
                model.set_voltage_protocol(param)
                # Check if the protocol gives nonsense simulations
                if not np.all(np.isfinite(model.simulated_currents)):
                    return np.inf

            ds = 100  # NOTE downsample rate for the output time points

            # Run simulations
            scores = []
            for model in self._model_list:
                # Run simulations
                ## Run simulation per parameter sample
                #sims = []
                #for i, p in enumerate(self._parameter_samples):
                #    sims.append(model.simulate(p))
                #sims = np.asarray(sims)
                # Broadcast method to run simulations for all parameter samples
                sims = model.simulate(self._parameter_samples, downsample=ds,
                        multi_input=True)

                if not np.all(np.isfinite(sims)):
                    return np.inf

                ds_time_points = len(model.times()[::ds])
                S_gsa = np.zeros((ds_time_points, self._n_model_parameters))
                ## Use SALib.analyze.sobol
                #for i, idx in enumerate(time_idx):
                #    # Compute Sobol sensitivity
                #    Si = sobol.analyze(self._salib_problem, sims[:, idx],
                #            calc_second_order=False)
                #    # Schenkendorf et al. 2018 (Eq. 10) uses only S1
                #    S_gsa[i, :] = Si['S1']
                # Use reimplemented sobol
                Si = sobol.analyze(self._salib_problem, sims, conf_level=None,
                        calc_second_order=False)

                S_gsa[:, :] = Si['S1'].T

                ds_t = model.times()[::100]  # time point for each S_gsa

                s1 = shannon_all_measure(S_gsa)
                s2 = shannon_segments_measure(S_gsa, ds_t, seg_t)
                s3 = parameter_dependency_measure(S_gsa)
                s4 = output_uncertainty_measure(sims.T)
                s = np.dot([s1, s2, s3, s4], self._weight)

                scores.append(s)

            # Take mean over all input models
            s = np.mean(scores)

            if np.isfinite(s):
                return s
            else:
                return np.inf
        except:
            return np.inf  # Just make sure it won't break the opt process.
