#
# PINTS and PyOED Model.
#
import os
import numpy as np
import pints
import pyoed
import myokit

from method.default import *

vhold = -80  # mV

parameters = [
    'ina.s', 'ical.s', 'ikr.s', 'iks.s', 'ito.s', 'inaca.s', 'ik1.s', 'inak.s'
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
        self.b = myokit.tools.Benchmarker()
    def exit(self):
        pass
    def update(self, progress):
        return self.b.time() < self.max_time


#
# Create ForwardModel
#

class CCModel(pints.ForwardModel, pyoed.ForwardModel):
    """
    # A current clamp (CC) model linking Myokit and Pints, PyOED ForwardModel.
    """

    def __init__(self, model_file, prepace=10, transform=None, dt=0.1,
                 norm=False, parameters=parameters, n_steps=None,
                 max_evaluation_time=60):
        """
        # model_file: mmt model file for myokit; main units: mV, ms, pA.
        # prepace: number of pre-pace before recording the simulated AP.
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
        self.dt = dt
        self._n_steps = n_steps
        #self._stimulate = stimulate
        #self._stim_seq = stim_seq
        self.transform = transform
        self.norm = norm  # if normalise
        self.max_evaluation_time = max_evaluation_time  # maximum time allowed
        self.parameters = parameters
        self._voltage_name = model_voltage[self._model_file_name]

        self.presimulation = myokit.Simulation(self._model)
        self.simulation = myokit.Simulation(self._model)

        # self.presimulation.set_tolerance(1e-8, 1e-10)
        self.presimulation.set_max_step_size(1e-1)  # ms
        # self.simulation.set_tolerance(1e-8, 1e-10)
        self.simulation.set_max_step_size(1e-1)  # ms

        ## Set up patch clamp
        #for ion_var, ion_conc in model_ion[self._model_file_name]:
        #    self._fix_concentration(self._model, ion_var, ion_conc)

        # Set stimulus default level
        try:
            stim_amp_var, stim_amp_val = model_stim_amp[self._model_file_name]
            self.presimulation.set_constant(stim_amp_var, stim_amp_val)
            self.simulation.set_constant(stim_amp_var, stim_amp_val)
        except:
            raise ValueError('Model stimulus do not exist in the given ' \
                             + 'model')

        # Add prepace stimulus
        self._stim_dur, stim_offset, cl, self._stim_amp = \
                model_stim_setup[self._model_file_name]
        self._prepace_cl = cl
        preprotocol = myokit.pacing.blocktrain(period=self._prepace_cl,
                                               duration=self._stim_dur, 
                                               offset=stim_offset,
                                               level=self._stim_amp)
        self.presimulation.set_protocol(preprotocol)
        del(preprotocol)

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
        self.use_init_state(False)
        self.set_init_state(self.original_state)
        print('Done')

    ############################################################################
    # Interface to PINTS and PyOED
    ############################################################################
    def n_parameters(self):
        return self.n_model_parameters()

    def n_model_parameters(self):
        return len(self.parameters)

    def n_design_variables(self):
        return self._n_steps * 2

    def design(self, variables):
        # Assume protocol p is
        # [stim_1_amp, holding_1_duration, stim_2_amp, ...]
        protocol = myokit.Protocol()
        duration = 0
        for i in range(len(variables) // 2):
            l = variables[2 * i]
            d = variables[2 * i + 1]
            # Add a stimulus at l size
            protocol.add_step(l * self._stim_amp, self._stim_dur)
            # Add holding (no stimulus) for d duration
            protocol.add_step(0, d)
            duration += self._stim_dur + d
        self.simulation.set_protocol(protocol)
        del(protocol)
        self._times = np.arange(0, duration, self.dt)

    def simulate(self, parameter, times=None, extra_log=[]):
        """
        Generate APs using the given scalings.
        """
        parameter = np.array(parameter)

        # Update model parameters
        if self.transform is not None:
            parameter = self.transform(parameter)
        # Simulate with modified model
        for i, name in enumerate(self._conductance):
            self.presimulation.set_constant(name,
                    parameter[i] * self.original[i])
            self.simulation.set_constant(name,
                    parameter[i] * self.original[i])

        # Set states
        if self._use_init_state:
            # Use a previously simulated states
            state_to_set = self._init_state
        else:
            self.presimulation.reset()
            # As myokit.org specified, in AP simulation mode, simulation.pre()
            # sorts the end of simulation state as the new default state, so
            # simulation.reset() only reset to the 'new' default state. It need
            # a manual reset of the state using simulation.set_state() to the
            # originally stored state.
            self.presimulation.set_state(self.original_state)
            try:
                # Pre-pace for some beats
                self.presimulation.pre(self._prepace * self._prepace_cl)
            except (myokit.SimulationError, myokit.SimulationCancelledError):
                return np.ones(self._times.shape) * float('inf')
            # Get the state from pre-pacing
            state_to_set = self.presimulation.state()
        self.simulation.reset()
        self.simulation.set_state(state_to_set)

        # Run simulation
        try:
            # Log some beats
            p = Timeout(self.max_evaluation_time)
            d = self.simulation.run(np.max(self._times)+1e-2,
                log_times=self._times,
                log=[self._voltage_name] + extra_log,
                progress=p,
                ).npview()
            del(p)
        except (myokit.SimulationError, myokit.SimulationCancelledError):
            return np.ones(self._times.shape) * float('inf')

        if self.norm:
            d[self._voltage_name] = self.normalise(d[self._voltage_name])

        if extra_log:
            return d
        else:
            return d[self._voltage_name]
    ############################################################################

    def set_init_state(self, states):
        self._init_state = states

    def use_init_state(self, use=True):
        self._use_init_state = use

    def _fix_concentration(self, model, variable, concentration):
        v = model.get(variable)
        if v.is_state():
            v.demote()
        v.set_rhs(concentration)
    
    def parameter(self):
        # return the name of the parameters
        return self.parameters

    def times(self):
        # Return the time series for the currently set stimulus protocol.
        return self._times

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
        # NOTE: OK to pre-compute the init states for VC when we change only
        #       the maximum conductance of the model!
        self.default_init_state = self.simulation1.state()
        self.simulation1.pre(100)
        self.set_init_state(self.simulation1.state())

        # Create simulation protocol
        self.simulation2 = myokit.Simulation(self._model)
        p = [self._vhold, 100] * 3 * len(self.current)
        self.design(p)
        # self.simulation2.set_tolerance(1e-6, 1e-8)
        # self.simulation2.set_max_step_size(1e-2)  # ms

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

    ############################################################################
    # Interface to PINTS and PyOED
    ############################################################################
    def n_parameters(self):
        return self.n_model_parameters()

    def n_model_parameters(self):
        return len(self.parameters)

    def n_design_variables(self):
        return self._n_steps * 2

    def design(self, variables):
        # self.set_voltage_protocol(variables)
        # Assume protocol p is
        # [step_1_voltage, step_1_duration, step_2_voltage, ...]
        # prt_mask: (numpy) mask function that remove part of the measurement;
        #           can be used as a capacitive filter, or to make the fitting
        #           harder
        protocol = myokit.Protocol()
        duration = 0
        p = variables
        for i in range(len(p) // 2):
            protocol.add_step(p[2 * i], p[2 * i + 1])
            duration += p[2 * i + 1]
        self.simulation2.set_protocol(protocol)
        del(protocol)
        self._times = np.arange(0, duration, self.dt)
        self.simulated_currents = self._simulate_protocol(self._times)

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
    ############################################################################

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
        self.simulation2.reset()
        self.simulation2.set_state(self.init_state)

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
        self.simulation2.reset()
        self.simulation2.set_state(self.init_state)

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
