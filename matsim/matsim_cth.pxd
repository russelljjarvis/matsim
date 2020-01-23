from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "bio.cpp":
    pass

# Decalre the class with cdef
cdef extern from "bio.h":
    cdef cppclass Conductance:
        void set_rate(double)
        double get_g()

    cdef cppclass ExponentialConductance(Conductance):
        ExponentialConductance(double, double, double) except +
        double get_g()
        double reversal

    cdef cppclass ShotNoiseConductance(Conductance):
        # ShotNoiseConductance() except +
        ShotNoiseConductance(double, double, double, double) except +
        void set_rate(double)
        double get_g()
        double reversal

    cdef cppclass OUConductance(Conductance):
        # ShotNoiseConductance() except +
        OUConductance(double, double, double, double) except +
        void set_rate(double)
        double get_g()
        double reversal

    cdef cppclass MATThresholds:
        MATThresholds() except +
        MATThresholds(double, double, double, double, double, double, bool) except +
        vector[double] get_spike_times()
        void reset_spike_times()
        double threshold
        bool resetting

    cdef cppclass Neuron:
        Neuron() except +
        Neuron(double, double, double, vector[MATThresholds*]) except +
        void append_conductance(Conductance)
        void timestep(double)
        vector[Conductance*] conductances
        vector[MATThresholds*] mats
        double voltage, time

    cdef cppclass HHNeuron:
        HHNeuron() except +
        HHNeuron(double) except +
        void append_conductance(Conductance)
        void timestep(double)
        vector[Conductance*] conductances
        double V, time
        double i_na, i_k, i_l, i_m
        double m, h, n, p

cdef extern from "simulation.cpp":
    vector[int] sr_experiment(Neuron neuron, double time_window, double dt,
                       vector[double] exc_intensities, vector[double] inh_intensities, int seed)


    vector[vector[double]] sr_experiment_spike_times(Neuron neuron, double time_window, double dt,
                       vector[double] exc_intensities, vector[double] inh_intensities, int seed)
