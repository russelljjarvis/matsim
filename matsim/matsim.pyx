# distutils: language = c++

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
from cython.operator cimport dereference as deref
from .matsim_cth cimport ExponentialConductance as CExponentialConductance
from .matsim_cth cimport ShotNoiseConductance as CShotNoiseConductance
from .matsim_cth cimport OUConductance as COUConductance
from .matsim_cth cimport Conductance as CConductance
from .matsim_cth cimport MATThresholds as CMATThresholds
from .matsim_cth cimport Neuron as CNeuron
from .matsim_cth cimport HHNeuron as CHHNeuron
from .matsim_cth cimport sr_experiment as _sr_experiment
from .matsim_cth cimport sr_experiment_spike_times as _sr_experiment_spike_times

import numpy as np
import pandas as pd

# Create a Cython extension type which holds a C++ instance
# as an attribute and create a bunch of forwarding methods
# Python extension type.

cdef class Conductance:
    cdef CConductance* conductance

    def set_rate(self, double rate):
        deref(self.conductance).set_rate(rate)

    @property
    def g(self):
        return deref(self.conductance).get_g()

    def set_g(self, g):
        deref(self.conductance).set_g(g)

cdef class ExponentialConductance(Conductance):
    cdef double g_peak, reversal, decay

    def __cinit__(self, double g_peak, double reversal, double decay):
        self.g_peak = g_peak
        self.reversal = reversal
        self.decay = decay
        self.conductance = new CExponentialConductance(g_peak, reversal, decay)

    def __dealloc__(self):
        del self.conductance

    def copy(self):
        new_conductance = ExponentialConductance(self.g_peak, self.reversal, self.decay)
        new_conductance.set_g(self.g)
        return new_conductance

cdef class ShotNoiseConductance(Conductance):
    cdef double rate, g_peak, reversal, decay
    # cdef CShotNoiseConductance* conductance  # Hold a C++ instance which we're wrapping

    def __cinit__(self, double rate, double g_peak, double reversal, double decay):
        self.rate = rate
        self.g_peak = g_peak
        self.reversal = reversal
        self.decay = decay
        self.conductance = new CShotNoiseConductance(rate, g_peak, reversal, decay)

    def __dealloc__(self):
        del self.conductance

    def copy(self):
        new_conductance = ShotNoiseConductance(self.rate, self.g_peak, self.reversal, self.decay)
        new_conductance.set_g(self.g)
        return new_conductance

    # def set_rate(self, double rate):
    #     deref(self.conductance).set_rate(rate)

cdef class OUConductance(Conductance):
    cdef double rate, g_peak, reversal, decay
    # cdef COUConductance* conductance  # Hold a C++ instance which we're wrapping

    def __cinit__(self, double rate, double g_peak, double reversal, double decay):
        self.rate = rate
        self.g_peak = g_peak
        self.reversal = reversal
        self.decay = decay
        self.conductance = new COUConductance(rate, g_peak, reversal, decay)

    def __dealloc__(self):
        del self.conductance

    def copy(self):
        new_conductance = OUConductance(self.rate, self.g_peak, self.reversal, self.decay)
        new_conductance.set_g(self.g)
        return new_conductance

    # @property
    # def g(self):
    #     return deref(self.conductance).g

    # def set_rate(self, double rate):
    #     deref(self.conductance).set_rate(rate)

cdef class MATThresholds:
    cdef CMATThresholds* mat  # Hold a C++ instance which we're wrapping
    cdef string name

    cdef double alpha1, alpha2, tau1, tau2, omega, refractory_period
    cdef bool resetting

    def __cinit__(self, double alpha1, double alpha2, double tau1, double tau2, double omega,
            double refractory_period, name, resetting=False):
        self.mat = new CMATThresholds(alpha1, alpha2, tau1, tau2, omega,
            refractory_period, resetting)
        self.name = <string> name.encode('utf-8')

        self.alpha1 = alpha1
        self.alpha2 = alpha2
        self.tau1   = tau1
        self.tau2   = tau2
        self.omega  = omega
        self.refractory_period = refractory_period
        self.resetting = resetting

    def __dealloc__(self):
        del self.mat

    def copy(self):
        new_mat = Neuron(self.alpha1, self.alpha2, self.tau1, self.tau2, self.omega,
            self.refractory_period, self.name, self.resetting)
        return new_mat

    @property
    def threshold(self):
        return deref(self.mat).threshold

    def get_spike_times(self):
        cdef vector[double] spike_times

        spike_times = deref(self.mat).get_spike_times()
        return np.array([ x for x in spike_times ])

    def reset_spike_times(self):
        deref(self.mat).reset_spike_times()


cdef class Neuron:
    cdef CNeuron neuron
    cdef vector[string] mat_names
    cdef double resting_potential, membrane_resistance, membrane_capacitance
    cdef object thresholds

    def __cinit__(self, double resting_potential, double membrane_resistance,
        double membrane_capacitance, mats):
        # cdef MATThresholds* c_mat
        cdef MATThresholds mat
        cdef vector[CMATThresholds*] mat_vec

        self.resting_potential = resting_potential
        self.membrane_resistance = membrane_resistance
        self.membrane_capacitance = membrane_capacitance
        self.thresholds = []

        for mat in mats:
            self.thresholds.append(mat)

            mat_vec.push_back(mat.mat)
            self.mat_names.push_back(mat.name)

        self.neuron = CNeuron(resting_potential, membrane_resistance, membrane_capacitance, mat_vec)
        # self.mats = mats

    def append_conductance(self, Conductance cond):
        self.neuron.conductances.push_back(cond.conductance)

    cpdef void timestep(self, double dt):
        self.neuron.timestep(dt)

    def copy(self):
        new_mats = [mat.copy() for mat in self.thresholds]
        new_neuron = Neuron(self.resting_potential, self.membrane_resistance, self.membrane_capacitance,
            )
        return new_neuron

    # Attribute access
    @property
    def voltage(self):
        return self.neuron.voltage

    @property
    def time(self):
        return self.neuron.time
    @time.setter
    def time(self, time):
        self.neuron.time = time


cdef class HHNeuron:
    cdef CHHNeuron neuron
    cdef double adaptation, VS, Ah
    cdef object conductances

    def __cinit__(self, adaptation=0.07, VS=-10, Ah=0.128, m=0, h=0, n=0, p=0, V=200):
        self.adaptation = adaptation
        self.VS = VS
        self.Ah = Ah
        self.neuron = CHHNeuron(adaptation, VS, Ah, m, h, n, p, V)
        self.conductances = []
        # self.mats = mats

    def append_conductance(self, Conductance cond):
        self.conductances.append(cond)
        self.neuron.conductances.push_back(cond.conductance)

    cpdef void timestep(self, double dt):
        self.neuron.timestep(dt)

    def copy(self):
        new_neuron = HHNeuron(self.adaptation, self.VS, self.Ah, *self.gate_vars, self.voltage)
        return new_neuron

    # Attribute access
    @property
    def conductance(self):
        return self.conductances

    @property
    def voltage(self):
        return self.neuron.V

    @property
    def i_na(self):
        return self.neuron.i_na

    @property
    def i_k(self):
        return self.neuron.i_k

    @property
    def i_l(self):
        return self.neuron.i_l

    @property
    def i_m(self):
        return self.neuron.i_m

    @property
    def gate_vars(self):
        return np.array([
            self.neuron.m, self.neuron.h, self.neuron.n, self.neuron.p
            ])

    @property
    def time(self):
        return self.neuron.time
    @time.setter
    def time(self, time):
        self.neuron.time = time

def sr_experiment(Neuron neuron, double time_window, double dt,
        intensities, intensity_freq_func, int seed):
    exc_intensities, inh_intensities = np.array([
        np.array(intensity_freq_func(i))
            for i in intensities]).T * dt

    cdef CNeuron c_neuron = neuron.neuron
    cdef vector[double] c_exc = exc_intensities
    cdef vector[double] c_inh = inh_intensities
    cdef vector[int] results

    mat_names = [name.decode("utf-8") for name in neuron.mat_names]

    results = _sr_experiment(c_neuron, time_window, dt, c_exc, c_inh, seed)
    result_array = np.array([x for x in results])

    return pd.DataFrame(result_array.reshape(-1, len(mat_names)), columns=mat_names, index=intensities).\
            groupby(level=0).agg(lambda x: list(x)).stack().swaplevel()

def sr_experiment(Neuron neuron, time_windows, dt, intensities, intensity_freq_func, seed):
    exc_intensities, inh_intensities = np.array([
        np.array(intensity_freq_func(i))
            for i in intensities]).T

    cdef CNeuron c_neuron = neuron.neuron
    cdef vector[double] c_exc = exc_intensities
    cdef vector[double] c_inh = inh_intensities
    cdef vector[vector[double]] results
    cdef double st

    mat_names = [name.decode("utf-8") for name in neuron.mat_names]

    results = _sr_experiment_spike_times(c_neuron, max(time_windows), dt, c_exc, c_inh, seed)
    result_python = [
        np.array([st for st in spike_times]) for spike_times in results
    ]
    
    result_dict = {}

    for tw in time_windows:
        result_dict[tw] = pd.DataFrame(
            np.array([ (x < tw).sum() for x in result_python ]).reshape(-1, len(mat_names)),
            columns=mat_names,
            index=intensities).groupby(level=0).agg(lambda x: list(x))

    return result_dict

def steady_spike_train(Neuron neuron, double time, double dt, exc, inh):
    mat_names = [name.decode("utf-8") for name in neuron.mat_names]
    spike_trains = {}
    cdef vector[double] spike_times
    cdef CMATThresholds* mat
    cdef CConductance *conductance

    conductance = neuron.neuron.conductances[0]
    deref(conductance).set_rate(exc)

    conductance = neuron.neuron.conductances[1]
    deref(conductance).set_rate(inh)

    cdef double tot_time = 0
    while tot_time < time:
        neuron.timestep(dt)
        tot_time += dt

    for i, name in enumerate(mat_names):
        mat = neuron.neuron.mats[i]
        spike_times = deref(mat).get_spike_times()
        deref(mat).reset_spike_times()
        spike_trains[name] = np.array([t for t in spike_times])

    return spike_trains
