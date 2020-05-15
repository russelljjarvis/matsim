#ifndef NEURON_H
#define NEURON_H
using namespace std;

class Conductance {
	double g, reversal;

	public:
		virtual void activate();
		virtual void update(double);
		virtual void set_rate(double);
		virtual double get_g();
		virtual void set_g(double);
		virtual double get_reversal();
};

class ExponentialConductance: public Conductance {
	double g_peak, decay, g, reversal;

	public:
		ExponentialConductance();
		ExponentialConductance(double, double, double);
		void update(double);
		double get_g();
		void set_g(double);
		double get_reversal();
		void activate();
};

class ShotNoiseConductance: public Conductance {
	double rate, g_peak, decay, g, reversal;

	public:
		ShotNoiseConductance();
		ShotNoiseConductance(double, double, double, double);
		void update(double);
		void set_rate(double);
		double get_g();
		void set_g(double);
		double get_reversal();
		void activate();
};

class OUConductance: public Conductance {
	double rate, g_peak, decay, g, reversal;
	double mean, sigma, D;
	double get_A(double);

	public:
		OUConductance();
		OUConductance(double, double, double, double);
		void update(double);
		void set_rate(double);
		double get_g();
		void set_g(double);
		double get_reversal();
		void activate();
};

class MATThresholds {
	double alpha1, alpha2, tau1, tau2, omega, t1, t2, refractory_period, past_spike_time;
	vector<double> spike_times;

	public:
		MATThresholds();
		MATThresholds(double, double, double, double, double, double, bool);
		void fire(double);
		void update(double);
		vector<double> get_spike_times();
		void reset_spike_times();
		double threshold;
		bool resetting;
};

class Neuron {
	double resting_potential, membrane_resistance, membrane_capacitance, time_constant, reset_potential;

	public:
		Neuron();
		Neuron(double, double, double, vector<MATThresholds*>, double);
		void append_conductance(Conductance*);
		void integrate_voltage(double);
		void timestep(double);
		vector<Conductance*> conductances;
		vector<MATThresholds*> mats;
		double voltage, time;
};

class MCNeuron {
	double resting_potential, membrane_resistance, leaky_conductance, membrane_capacitance, time_constant, reset_potential;
	double coupling_conductance;

	public:
		MCNeuron();
		MCNeuron(double, double, double, vector<MATThresholds*>, double, double);
		void append_conductance(Conductance*, int);
		void integrate_voltage(double);
		void timestep(double);
		vector<Conductance*> conductances_soma, conductances_dendrite;
		vector<MATThresholds*> mats;
		double voltageSoma, voltageDendrite, time;
};

class HHNeuron {
	double g_l, E_l, c_m, E_na, g_na, E_k, g_k, g_m, VS, Ah;

	public:
		HHNeuron();
		HHNeuron(double, double, double, double, double, double, double, double);
		void append_conductance(Conductance*);
		void integrate_voltage(double);
		void timestep(double);
		vector<Conductance*> conductances;
		double V, time;
		double i_na, i_k, i_l, i_m;
		double m, h, n, p;
};

#endif
