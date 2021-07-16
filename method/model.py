#
# PINTS Model.
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals
import os
import numpy as np
import pints
import pyoed
import myokit

from default import *

vhold = -80  # mV
DT = 0.2

#
# Defining AP simulation
#

def simulate_aps(scales, model_file, beats=2, cl=1000, prepace=100,
        stimulate=True):
    """
    # A simple method to generate APs using the given scalings.
    """
    # Load model
    model = myokit.load_model(model_file)

    # Apply scalings
    for var in scales:
        scale = scales[var]
        v = model.get(var)
        v.set_rhs(myokit.Multiply(myokit.Number(scale), v.rhs()))

    # Simulate with modified model
    sim = myokit.Simulation(model)

    # Add stimulus
    if stimulate:
        protocol = myokit.pacing.blocktrain(period=cl, duration=1,
                offset=50)
        sim.set_protocol(protocol)
    
    # Pre-pace for some beats
    sim.pre(prepace * cl)

    # Log some beats and return
    log = ['engine.time', 'membrane.V']
    log = ['environment.time', 'membrane.V']
    log = ['membrane.V']
    return sim.run(beats * cl, log=log, log_interval=DT).npview()


#
# Time out handler
#
class Timeout(myokit.ProgressReporter):
    """
    A :class:`myokit.ProgressReporter` that halts the simulation after
    ``max_time`` seconds.
    """
    def __init__(self, max_time):
        self.max_time = float(max_time)
    def enter(self, msg=None):
        self.b = myokit.Benchmarker()
    def exit(self):
        pass
    def update(self, progress):
        return self.b.time() < self.max_time


parameters = [
        'ina.s', 'ical.s', 'ikr.s', 'iks.s', 'ito.s', 'inaca.s', 'ik1.s',
        'inak.s', #'if.s',
        ]
parameter_names = [
        r'$s_{Na}$',
        r'$s_{CaL}$',
        r'$s_{Kr}$',
        r'$s_{Ks}$',
        r'$s_{to}$',
        r'$s_{NaCa}$',
        r'$s_{K1}$',
        r'$s_{NaK}$'
        ]


#
# Create ForwardModel
#

class CCModel(pints.ForwardModel):
    """
    # A current clamp (CC) model linking Myokit and Pints ForwardModel.
    """

    def __init__(self, model_file, prepace=10, stimulate=True, stim_seq=None,
            transform=None, max_evaluation_time=5, norm=False):
        """
        # model_file: mmt model file for myokit; main units: mV, ms, pA.
        # prepace: number of pre-pace before recording the simulated AP.
        # stimulate: bool, if True, apply stimulus.
        # stim_seq: array-like, a sequence of stimulus in
        #           [(stim_level1, duration1), (stim_level2, duration2)...]
        #           e.g. [(0, 50), (1, 5), (0, 945)]
        # transform: transform search space parameters to model parameters.
        # max_evaluation_time: maximum time (in second) allowed for one
        #                      simulate() call.
        # norm: bool, if True, normalise output.
        """
        self._model = myokit.load_model(model_file)
        self._model_file = model_file
        self._model_file_name = os.path.basename(model_file)
        print('Initialising model %s...' % self._model_file_name)
        self._prepace = prepace
        self._stimulate = stimulate
        self._stim_seq = stim_seq
        self.transform = transform
        self.presimulation = myokit.Simulation(self._model)
        self.simulation = myokit.Simulation(self._model)
        self.parameters = parameters

        # self.presimulation.set_tolerance(1e-8, 1e-10)
        # self.presimulation.set_max_step_size(1e-2)  # ms
        # self.simulation.set_tolerance(1e-8, 1e-10)
        # self.simulation.set_max_step_size(1e-2)  # ms

        # Set stimulus default level
        try:
            stim_amp_var, stim_amp_val = model_stim_amp[self._model_file_name]
            self.presimulation.set_constant(stim_amp_var, stim_amp_val)
            self.simulation.set_constant(stim_amp_var, stim_amp_val)
        except:
            raise ValueError('Model stimulus do not exist in the given ' \
                    + 'model')

        # Add prepace stimulus
        stim_dur, stim_offset, cl, stim_amp = \
                model_stim_setup[self._model_file_name]
        self._prepace_cl = cl
        if '-old' in self._model_file_name:
            self._prepace_cl = self._prepace_cl * 1e-3
            stim_dur = stim_dur * 1e-3
            stim_offset = stim_offset * 1e-3
        preprotocol = myokit.pacing.blocktrain(period=self._prepace_cl,
                                               duration=stim_dur, 
                                               offset=stim_offset,
                                               level=stim_amp)
        self.presimulation.set_protocol(preprotocol)
        del(preprotocol)

        # Add stimulus
        if self._stimulate:
            if stim_seq is not None:
                protocol = myokit.Protocol()
                for l, t in self._stim_seq:
                    if l > 0:
                        protocol.add_step(l * stim_amp, stim_dur)
                    else:
                        protocol.add_step(l, t)
                self.simulation.set_protocol(protocol)
            else:
                protocol = myokit.pacing.blocktrain(period=cl,
                                                    duration=stim_dur, 
                                                    offset=stim_offset,
                                                    level=stim_amp)
                self.simulation.set_protocol(protocol)
            del(protocol)

        # Create a order-matched conductance list
        try:
            self._conductance = []
            p2g = model_conductance[self._model_file_name]
            for p in self.parameters:
                self._conductance.append(p2g[p])
            assert(len(self._conductance) == len(self.parameters))
            del(p, p2g)
        except:
            raise ValueError('Model conductances do not match parameters')

        # Get original parameters
        try:
            self.original = []
            for name in self._conductance:
                v = self._model.get(name)
                self.original.append(np.float(v.rhs()))
            assert(len(self.original) == len(self.parameters))
        except:
            raise ValueError('Model conductances do not exist in the given ' \
                    + 'model')

        # Store model original state -- can be something else later!
        self.original_state = self.simulation.state()
        # if normalise
        self.norm = norm
        # maximum time allowed
        self.max_evaluation_time = max_evaluation_time
        print('Done')
    
    def n_parameters(self):
        return len(self.parameters)

    def simulate(self, parameter, times, extra_log=[]):
        """
        Generate APs using the given scalings.
        """
        parameter = np.array(parameter)

        if '-old' in self._model_file_name:
            times = times * 1e-3

        # Update model parameters
        if self.transform is not None:
            parameter = self.transform(parameter)
        # Simulate with modified model
        for i, name in enumerate(self._conductance):
            self.presimulation.set_constant(name,
                    parameter[i] * self.original[i])
            self.simulation.set_constant(name,
                    parameter[i] * self.original[i])

        # Run
        self.presimulation.reset()
        self.simulation.reset()
        # As myokit.org specified, in AP simulation mode, simulation.pre()
        # sorts the end of simulation state as the new default state, so
        # simulation.reset() only reset to the 'new' default state. It need
        # a manual reset of the state using simulation.set_state() to the 
        # originally stored state.
        self.presimulation.set_state(self.original_state)
        try:
            # Pre-pace for some beats
            self.presimulation.pre(self._prepace * self._prepace_cl)
            self.simulation.set_state(self.presimulation.state())
            # Log some beats
            p = Timeout(self.max_evaluation_time)
            d = self.simulation.run(np.max(times)+0.02, 
                log_times=times,
                log=['membrane.V'] + extra_log,
                progress=p,
                ).npview()
            del(p)
        except (myokit.SimulationError, myokit.SimulationCancelledError):
            return np.ones(times.shape) * float('inf')

        if self.norm:
            d['membrane.V'] = self.normalise(d['membrane.V'])


        if '-old' in self._model_file_name:
            d['membrane.V'] = d['membrane.V'] * 1e3

        if extra_log:
            return d
        else:
            return d['membrane.V']

    def _fix_concentration(self, model, variable, concentration):
        v = model.get(variable)
        if v.is_state():
            v.demote()
        v.set_rhs(concentration)
    
    def parameter(self):
        # return the name of the parameters
        return self.parameters

    def name(self):
        # name
        return self._name
    
    def set_name(self, name):
        # set name
        self._name = name

    def normalise(self, v):
        # Do whatever fancy normalisation here
        # For example to mimic optical mapping data
        method = 1
        if method == 1:
            # 5, 95 percentiles
            minimum = np.percentile(v, 5)
            maximum = np.percentile(v, 95)
            return (v - minimum) / (maximum - minimum)
        elif method == 2:
            # RMSD normalisation
            return v / np.sqrt(np.mean(v ** 2))
        elif method == 3:
            # Use minimisation to fit the two curve
            # Actually it should not be here but in
            # the error function...
            raise NotImplementedError


class VCModel(pints.ForwardModel, pyoed.ForwardModel):
    """
    # A voltage clamp (VC) model linking Myokit and Pints, PyOED ForwardModel.
    """

    def __init__(self, model_file, transform=None, dt=0.1,
            parameters=parameters, n_steps=None, max_evaluation_time=60):
        """
        # model_file: mmt model file for myokit; main units: mV, ms, pA.
        # transform: transform search space parameters to model parameters.
        # dt: sample duration in ms.
        # max_evaluation_time: maximum time (in second) allowed for one
        #                      simulate() call.
        """
        self._model = myokit.load_model(model_file)
        self._model_file = model_file
        self._model_file_name = os.path.basename(model_file)
        print('Initialising model %s...' % self._model_file_name)
        self.transform = transform
        self.parameters = parameters
        self.dt = dt
        self._vhold = vhold
        self._n_steps = n_steps
        self._voltage_name = model_voltage[self._model_file_name]

        # maximum time allowed
        self.max_evaluation_time = max_evaluation_time

        # Get current names of output
        self.current = []
        m_cur = model_current[self._model_file_name]
        for name in self.parameters:
            self.current.append(m_cur[name])
        # Set up voltage clamp
        for ion_var, ion_conc in model_ion[self._model_file_name]:
            self._fix_concentration(self._model, ion_var, ion_conc)
        # Detach voltage for voltage clamp(?)
        model_v = self._model.get(self._voltage_name)
        model_v.demote()
        model_v.set_rhs(self._vhold)
        self._model.get(model_pace[self._model_file_name]).set_binding(None)
        model_v.set_binding('pace')

        # Create pre-pacing protocol
        protocol = myokit.pacing.constant(self._vhold)
        # Create pre-pacing simulation
        self.simulation1 = myokit.Simulation(self._model, protocol)

        # Init states
        self.default_init_state = self.simulation1.state()
        self.init_state = self.default_init_state

        # Create simulation protocol
        self.simulation2 = myokit.Simulation(self._model)
        p = [self._vhold, 100] * 3 * len(self.current)
        self.set_voltage_protocol(p)
        self.simulation2.set_tolerance(1e-6, 1e-8)
        self.simulation2.set_max_step_size(1e-2)  # ms

        # Create a order-matched conductance list
        try:
            self._conductance = []
            p2g = model_conductance[self._model_file_name]
            for p in self.parameters:
                self._conductance.append(p2g[p])
            assert(len(self._conductance) == len(self.parameters))
            del(p, p2g)
        except:
            raise ValueError('Model conductances do not match parameters')

        # Get original parameters
        try:
            self.original = []
            for name in self._conductance:
                v = self._model.get(name)
                self.original.append(np.float(v.rhs()))
            assert(len(self.original) == len(self.parameters))
        except:
            raise ValueError('Model conductances do not exist in the given ' \
                    + 'model')
        print('Done')

    ###########################################################################
    # Interface to PINTS and PyOED
    ###########################################################################
    def n_parameters(self):
        return self.n_model_parameters()

    def n_model_parameters(self):
        return len(self.parameters)

    def n_design_variables(self):
        return self._n_steps * 2

    def design(self, variables):
        self.set_voltage_protocol(variables)

    def simulate(self, parameter, times=None, downsample=None,
                 multi_input=False):
        """
        Generate current of voltage clamp using the given scalings.

        Pre-simulated each current when setting the voltage clamp protocol.
        Time is a dummy argument for PINTS only.

        downsample: (int) downsample rate for the output time series.
        multi_input: (bool) if True, the input `parameter` is a list of
                     parameters to be simulated, and return a list of output
                     matching the order of the input parameters.
        """
        parameter = np.array(parameter)
        # Update model parameters
        if self.transform is not None:
            parameter = self.transform(parameter)
        if multi_input:
            # parameter = [set1_parameters, set2_parameters, ...]
            # Broadcast `parameter` to
            # (n_parameter_sets, n_currents [i.e. n_param_per_set],
            #  1 [broadcast with n_t])
            n_ps, n_c = parameter.shape
            parameter = parameter.reshape((n_ps,)+(n_c,)+(1,))
            # Simulation
            sim = self.simulated_currents[:, :]
            if downsample is not None:
                # downsample time points
                sim = sim[:, ::downsample]
            # (n_ps, n_c, 1) * (n_c, n_t) -> (n_ps, n_c, n_t)
            # sum across axis=1 -> (n_ps, n_t)
            return np.sum(parameter * sim, axis=1)

        else:
            # parameter = set1_parameters
            p = np.asarray(parameter).reshape(-1, 1)  # turn to column matrix
            # Simulation
            sim = self.simulated_currents[:, :]
            if downsample is not None:
                # downsample time points
                sim = sim[:, ::downsample]
            # (n_c, 1) * (n_c, n_t) -> (n_c, n_t)
            # sum across axis=0 -> (n_t,)
            return np.sum(p * sim, axis=0)
    ###########################################################################

    def set_init_state(self, v):
        self.init_state = v

    def current_state(self):
        return self.simulation2.state()

    def _fix_concentration(self, model, variable, concentration):
        v = model.get(variable)
        if v.is_state():
            v.demote()
        v.set_rhs(concentration)

    def set_parameters(self, parameter_dict):
        for name in parameter_dict:
            self.simulation1.set_constant(name, parameter_dict[name])
            self.simulation2.set_constant(name, parameter_dict[name])

    def set_voltage_protocol(self, p, prt_mask=None):
        # Assume protocol p is
        # [step_1_voltage, step_1_duration, step_2_voltage, ...]
        # prt_mask: (numpy) mask function that remove part of the measurement;
        #           can be used as a capacitive filter, or to make the fitting
        #           harder
        protocol = myokit.Protocol()
        duration = 0
        for i in range(len(p) // 2):
            protocol.add_step(p[2 * i], p[2 * i + 1])
            duration += p[2 * i + 1]
        self.simulation2.set_protocol(protocol)
        del(protocol)
        self.prt_mask = prt_mask
        self._times = np.arange(0, duration, self.dt)
        if self.prt_mask is not None:
            self._times = self._times[self.prt_mask]
        self.simulated_currents = self._simulate_protocol(self._times)

    def set_fixed_form_voltage_protocol(self, v, t, prt_mask=None):
        # v, t: voltage, time to be set in ms, mV
        # prt_mask: (numpy) mask function that remove part of the measurement;
        #           can be used as a capacitive filter, or to make the fitting
        #           harder
        self.simulation2.set_fixed_form_protocol(
            t, v  # ms, mV
        )
        self.prt_mask = prt_mask
        #self._times = t
        duration = t[-1]
        self._times = np.arange(0, duration, self.dt)
        if self.prt_mask is not None:
            self._times = self._times[self.prt_mask]
        self.simulated_currents = self._simulate_protocol(self._times)

    def times(self):
        # Return the time series for the currently set voltage protocol.
        return self._times

    def _simulate_protocol(self, times):
        self.simulation1.reset()
        self.simulation2.reset()
        self.simulation1.set_state(self.init_state)
        self.simulation1.pre(100)
        self.simulation2.set_state(self.simulation1.state())

        try:
            p = Timeout(self.max_evaluation_time)
            d = self.simulation2.run(np.max(times) + 0.02,
                log_times=times,
                log=self.current,
                progress=p,
                ).npview()
            del(p)
            # Convert dict to 2d array
            o = []
            for c in self.current:
                o.append(np.asarray(d[c]))
            o = np.asarray(o)
        except (myokit.SimulationError, myokit.SimulationCancelledError):
            o = np.ones((len(self.parameters), len(times))) * float('inf')
        return o

    def current_list(self):
        return self.current

    def voltage(self, times):
        self.simulation1.reset()
        self.simulation2.reset()
        self.simulation1.set_state(self.init_state)
        self.simulation1.pre(100)
        self.simulation2.set_state(self.simulation1.state())

        try:
            p = Timeout(self.max_evaluation_time)
            d = self.simulation2.run(np.max(times) + 0.02,
                log_times=times,
                log=[self._voltage_name],
                progress=p,
                ).npview()
            del(p)
        except (myokit.SimulationError, myokit.SimulationCancelledError):
            d = {self._voltage_name: np.ones(times.shape) * float('inf')}
        return d[self._voltage_name]


class BootstrapVCModel(pints.ForwardModel):
    """
    # Bootstrapping voltage clamp (VC) models, linking Myokit and
    # Pints ForwardModel.
    """

    def __init__(self, model_files, n_samples=100, transform=None, dt=0.1,
            max_evaluation_time=60):
        """
        # model_files: a list of mmt model files for myokit;
        #              main units: mV, ms, pA.
        # n_samples: number of bootstrap resampling.
        # transform: transform search space parameters to model parameters.
        # dt: sample duration in ms.
        # max_evaluation_time: maximum time (in second) allowed for one
        #                      simulate() call.
        """
        self._models = [myokit.load_model(m) for m in model_files]
        self._model_files = model_files
        self._model_file_names = [os.path.basename(m) for m in model_files]
        print('Initialising model %s...' % self._model_file_names)
        self.transform = transform
        self.parameters = parameters

        self._n_models = len(self._model_files)
        self._n_samples = int(n_samples)
        self._n_parameters = len(self.parameters)

        # Common dt and vhold, in ms and mV
        self.dt = dt
        self._vhold = vhold
        # A list of dt and vhold for each model (in case different unit!)
        self._dts = [self.dt] * self._n_models
        self._vholds = [self._vhold] * self._n_models

        # Samples for bootstrapping models
        self._sample_indices = self._resample(self._n_samples,
                self._n_parameters, self._n_models)

        # maximum time allowed
        self.max_evaluation_time = max_evaluation_time

        self.currents = []
        self.simulation1s = []
        self.simulation2s = []
        self.default_init_states =[]
        self.init_states = []
        self._conductances = []
        self.originals = []
        for i, m in enumerate(self._models):
            model_file_name = self._model_file_names[i]

            if '-old' in model_file_name:
                self._dts[i] = self.dt * 1e-3
                self._vholds[i] = self._vhold * 1e-3

            # Get current names of output
            c = []
            m_cur = model_current[model_file_name]
            for name in self.parameters:
                c.append(m_cur[name])
            # Set up voltage clamp
            for ion_var, ion_conc in model_ion[model_file_name]:
                self._fix_concentration(m, ion_var, ion_conc)
            # Detach voltage for voltage clamp(?)
            model_v = m.get('membrane.V')
            model_v.demote()
            model_v.set_rhs(self._vholds[i])
            m.get('engine.pace').set_binding(None)
            model_v.set_binding('pace')

            # Create pre-pacing protocol
            protocol = myokit.pacing.constant(self._vholds[i])
            # Create pre-pacing simulation
            s1 = myokit.Simulation(m, protocol)

            # Init states
            self.default_init_states.append(s1.state())
            self.init_states.append(s1.state())

            # Create simulation protocol
            s2 = myokit.Simulation(m)
            s2.set_tolerance(1e-6, 1e-8)
            s2.set_max_step_size(1e-2)  # ms

            self.currents.append(c)
            self.simulation1s.append(s1)
            self.simulation2s.append(s2)

            # Create a order-matched conductance list
            try:
                conductance = []
                p2g = model_conductance[model_file_name]
                for p in self.parameters:
                    conductance.append(p2g[p])
                assert(len(conductance) == self._n_parameters)
                self._conductances.append(conductance)
                del(p, p2g)
            except:
                raise ValueError('Model conductances do not match parameters')

            # Get original parameters
            try:
                original = []
                for name in conductance:
                    v = m.get(name)
                    original.append(np.float(v.rhs()))
                assert(len(original) == self._n_parameters)
                self.originals.append(original)
            except:
                raise ValueError('Model conductances do not exist in the ' \
                        + 'given model')

        # Set init protocol
        p = [self._vhold, 100] * 3 * len(c)
        self.set_voltage_protocol(p)

        print('Done')

    def n_parameters(self):
        return self._n_parameters

    def _resample(self, n_s, n_c, n_m):
        # Draw sample indices for the bootstrapping method.
        # n_s: number of samples.
        # n_c: number of model currents.
        # n_m: number of models.
        return np.random.randint(low=0, high=n_m, size=(n_s, n_c))

    def sample_indices(self):
        # Return the sample indices used in this bootstrapping.
        return self._sample_indices

    def set_init_state(self, v):
        # Set init states externally.
        # v: a list of initial states for each model
        self.init_states = v.copy()

    def current_state(self):
        # Return the a list of current states for each model
        return [s2.state() for s2 in self.simulation2s]

    def _fix_concentration(self, model, variable, concentration):
        v = model.get(variable)
        if v.is_state():
            v.demote()
        v.set_rhs(concentration)

    def set_voltage_protocol(self, prt, prt_mask=None):
        # Assume protocol `prt` is
        # [step_1_voltage, step_1_duration, step_2_voltage, ...] in ms, mV
        # prt_mask: (numpy) mask function that remove part of the measurement;
        #           can be used as a capacitive filter, or to make the fitting
        #           harder
        self.prt_mask = prt_mask
        self._times = np.arange(0, np.sum(np.asarray(prt)[1::2]), self.dt)
        if self.prt_mask is not None:
            self._times = self._times[self.prt_mask]
        self.simulated_currents = []
        for i_m in range(self._n_models):
            protocol = myokit.Protocol()
            if '-old' in self._model_file_names[i_m]:
                p = np.array(prt) * 1e-3
            else:
                p = np.array(prt)
            duration = 0
            for i in range(len(p) // 2):
                protocol.add_step(p[2 * i], p[2 * i + 1])
                duration += p[2 * i + 1]
            self.simulation2s[i_m].set_protocol(protocol)
            del(protocol)
            times = np.arange(0, duration, self._dts[i_m])
            if self.prt_mask is not None:
                times = times[self.prt_mask]
            self.simulated_currents.append(self._simulate_protocol(i_m, times))
        self.simulated_currents = np.asarray(self.simulated_currents)

    def set_fixed_form_voltage_protocol(self, v, t, prt_mask=None):
        # v, t: voltage, time to be set in ms, mV
        # prt_mask: (numpy) mask function that remove part of the measurement;
        #           can be used as a capacitive filter, or to make the fitting
        #           harder
        self.prt_mask = prt_mask
        self._times = t
        if self.prt_mask is not None:
            self._times = self._times[self.prt_mask]
        self.simulated_currents = []
        for i in range(self._n_models):
            if '-old' in self._model_file_names[i]:
                vi = v * 1e-3
                ti = t * 1e-3
            else:
                vi, ti = v, t
            self.simulation2s[i].set_fixed_form_protocol(
                ti, vi  # ms, mV
            )
            if self.prt_mask is not None:
                ti = ti[self.prt_mask]
            self.simulated_currents.append(self._simulate_protocol(i, ti))
        self.simulated_currents = np.asarray(self.simulated_currents)

    def times(self):
        # Return the time series for the currently set voltage protocol (ms).
        return self._times

    def _simulate_protocol(self, i_m, times):
        # i_m: the ith model
        # times: simulation time (in the right unit!)
        s1 = self.simulation1s[i_m]
        s2 = self.simulation2s[i_m]

        s1.reset()
        s2.reset()
        s1.set_state(self.init_states[i_m])
        s1.pre(100)
        s2.set_state(s1.state())

        try:
            p = Timeout(self.max_evaluation_time)
            d = s2.run(np.max(times) + 0.02,
                log_times=times,
                log=self.currents[i_m],
                progress=p,
                ).npview()
            del(p)
            # Convert dict to 2d array
            o = []
            for c in self.currents[i_m]:
                o.append(np.asarray(d[c]))
        except (myokit.SimulationError, myokit.SimulationCancelledError):
            o = np.ones((self._n_parameters, len(times))) * float('inf')
        return o

    def current_list(self):
        # Return a list of model's list of current name
        return self.currents

    def simulate(self, parameter, times=None, downsample=None,
                 multi_input=False, sample_index=None):
        """
        Generate current of voltage clamp using the given scalings.

        Pre-simulated each current when setting the voltage clamp protocol.
        Time is a dummy argument for PINTS only.

        downsample: (int) downsample rate for the output time series.
        multi_input: (bool) if True, the input `parameter` is a list of
                     parameters to be simulated, and return a list of output
                     matching the order of the input parameters.
        sample_index: (int) index of the 'sample model': `_sample_indices`.
        """
        parameter = np.array(parameter)
        # Update model parameters
        if self.transform is not None:
            parameter = self.transform(parameter)

        if multi_input and (sample_index is not None):
            # parameter = [set1_parameters, set2_parameters, ...]
            # Broadcast `parameter` to
            # (n_parameter_sets, n_currents [i.e. n_param_per_set],
            #  1 [broadcast with n_t])
            n_ps, n_c = parameter.shape
            parameter = parameter.reshape((n_ps,)+(n_c,)+(1,))
            r = np.arange(self._n_parameters)
            s = self._sample_indices[sample_index]
            sim = self.simulated_currents[s, r, :]
            if downsample is not None:
                # downsample time points
                sim = sim[:, ::downsample]
            # (n_ps, n_c, 1) * (n_c, n_t) -> (n_ps, n_c, n_t)
            # sum across axis=1 -> (n_ps, n_t)
            return np.sum(parameter * sim, axis=1)

        elif multi_input and (sample_index is None):
            # parameter = [set1_parameters, set2_parameters, ...]
            # Broadcast `parameter` to
            # (n_parameter_sets, 1 [broadcast with n_bs_samples],
            #  n_currents [i.e. n_param_per_set], 1 [broadcast with n_t])
            n_ps, n_c = parameter.shape
            parameter = parameter.reshape((n_ps,)+(1,)+(n_c,)+(1,))
            # Resample (mix) model currents with bootstrap indices
            r = np.arange(self._n_parameters)  # TODO could move this out
            sim = self.simulated_currents[self._sample_indices, r, :]
            if downsample is not None:
                sim = sim[:, :, ::downsample]
            # (n_ps, 1, n_c, 1) * (n_s, n_c, n_t) -> (n_ps, n_s, n_c, n_t)
            # sum across axis=2 -> (n_ps, n_s, n_t)
            sim = np.sum(parameter * sim, axis=2)
            return sim

        else:
            # parameter = set1_parameters
            p = np.asarray(parameter).reshape(-1, 1)  # turn to column matrix
            r = np.arange(self._n_parameters)
            out = []
            for i, s in enumerate(self._sample_indices):
                sim = self.simulated_currents[s, r, :]  # Get model and current
                if downsample is not None:
                    # downsample time points
                    sim = sim[:, ::downsample]
                out.append(np.sum(p * sim, axis=0))
            return np.asarray(out)

    def voltage(self, times):
        s1 = self.simulation1s[0]
        s2 = self.simulation2s[0]

        if '-old' in self._model_file_names[0]:
            times = times * 1e-3

        s1.reset()
        s2.reset()
        s1.set_state(self.init_states[0])
        s1.pre(100)
        s2.set_state(s1.state())

        try:
            p = Timeout(self.max_evaluation_time)
            d = s2.run(np.max(times) + 0.02,
                log_times=times,
                log=['membrane.V'],
                progress=p,
                ).npview()
            del(p)
        except (myokit.SimulationError, myokit.SimulationCancelledError):
            d = {'membrane.V': np.ones(times.shape) * float('inf')}

        if '-old' in self._model_file_names[0]:
            d['membrane.V'] = d['membrane.V'] * 1e-3

        return d['membrane.V']


class BMAVCModel(BootstrapVCModel):
    """
    # Bayesian model averaging (BMA) voltage clamp (VC) models, linking Myokit
    # and Pints ForwardModel.
    """

    def __init__(self, model_files, prior=None, transform=None, dt=0.1,
            max_evaluation_time=60):
        super(BMAVCModel, self).__init__(model_files, n_samples=1,
                transform=transform, dt=dt,
                max_evaluation_time=max_evaluation_time)
        del(self._n_samples, self._sample_indices)

        if prior is None:
            self._prior = np.ones(self._n_models) / np.float(self._n_models)
        else:
            self._prior = np.array(prior)
        self._prior = self._prior.reshape((-1, 1, 1))  # broadcast

    def simulate(self, parameter, times=None, downsample=None,
                 multi_input=False):
        """
        Generate current of voltage clamp using the given scalings.

        Pre-simulated each current when setting the voltage clamp protocol.
        Time is a dummy argument for PINTS only.

        downsample: (int) downsample rate for the output time series.
        multi_input: (bool) if True, the input `parameter` is a list of
                     parameters to be simulated, and return a list of output
                     matching the order of the input parameters.
        """
        parameter = np.array(parameter)
        # Update model parameters
        if self.transform is not None:
            parameter = self.transform(parameter)

        if multi_input:
            # parameter = [set1_parameters, set2_parameters, ...]
            # Broadcast `parameter` to
            # (n_parameter_sets, n_currents [i.e. n_param_per_set],
            #  1 [broadcast with n_t])
            n_ps, n_c = parameter.shape
            parameter = parameter.reshape((n_ps,)+(n_c,)+(1,))
            # Simulation and prior
            sim = self._prior * self.simulated_currents[:, :, :]
            sim = np.sum(sim, axis=0)
            if downsample is not None:
                # downsample time points
                sim = sim[:, ::downsample]
            # (n_ps, n_c, 1) * (n_c, n_t) -> (n_ps, n_c, n_t)
            # sum across axis=1 -> (n_ps, n_t)
            return np.sum(parameter * sim, axis=1)

        else:
            # parameter = set1_parameters
            p = np.asarray(parameter).reshape(-1, 1)  # turn to column matrix
            # Simulation and prior
            sim = self._prior * self.simulated_currents[:, :, :]
            sim = np.sum(sim, axis=0)
            if downsample is not None:
                # downsample time points
                sim = sim[:, ::downsample]
            # (n_c, 1) * (n_c, n_t) -> (n_c, n_t)
            # sum across axis=0 -> (n_t,)
            return np.sum(p * sim, axis=0)


class BMIAVCModel(BootstrapVCModel):
    """
    # Bayesian model 'independent' averaging (BMIA) voltage clamp (VC) models,
    # where (submodels) currents from each model (AP model) are treated
    # independently during fitting, linking Myokit and Pints ForwardModel.
    """

    def __init__(self, model_files, transform=None, dt=0.1,
            max_evaluation_time=60):
        super(BMIAVCModel, self).__init__(model_files, n_samples=1,
                transform=transform, dt=dt,
                max_evaluation_time=max_evaluation_time)
        del(self._n_samples, self._sample_indices)
        self._n_parameters *= self._n_models

    def simulate(self, parameter, times=None, downsample=None,
                 multi_input=False):
        """
        Generate current of voltage clamp using the given scalings.

        Pre-simulated each current when setting the voltage clamp protocol.
        Time is a dummy argument for PINTS only.

        downsample: (int) downsample rate for the output time series.
        multi_input: (bool) if True, the input `parameter` is a list of
                     parameters to be simulated, and return a list of output
                     matching the order of the input parameters.
        """
        parameter = np.array(parameter)
        # Update model parameters
        if self.transform is not None:
            parameter = self.transform(parameter)

        parameter = parameter.reshape(self._n_models, -1, 1)
        sim = self.simulated_currents[:, :, :]
        if downsample is not None:
            # downsample time points
            sim = sim[:, ::downsample]
        # (n_m, n_c, n_t) sum across axis=(0, 1) -> (n_t,)
        return np.sum(parameter * sim, axis=(0, 1))
