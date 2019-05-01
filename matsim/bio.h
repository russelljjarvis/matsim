#ifndef NEURON_H
#define NEURON_H
using namespace std;

class Conductance {
	double g, reversal;

	public:
		virtual void update(double);
		virtual void set_rate(double);
		virtual double get_g();
		virtual double get_reversal();
};

class ShotNoiseConductance: public Conductance {
	double rate, g_peak, decay, g, reversal;

	public:
		ShotNoiseConductance();
		ShotNoiseConductance(double, double, double, double);
		void update(double);
		void set_rate(double);
		double get_g();
		double get_reversal();
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
		double get_reversal();
};

class MATThresholds {
	double alpha1, alpha2, tau1, tau2, omega, t1, t2, refractory_period, past_spike_time;
	vector<double> spike_times;

	public:
		MATThresholds();
		MATThresholds(double, double, double, double, double, double);
		void fire(double);
		void update(double);
		vector<double> get_spike_times();
		void reset_spike_times();
		double threshold;
};

class Neuron {
	double resting_potential, membrane_resistance, membrane_capacitance, time_constant;

	public:
		Neuron();
		Neuron(double, double, double, vector<MATThresholds*>);
		void append_conductance(Conductance*);
		void integrate_voltage(double);
		void timestep(double);
		vector<Conductance*> conductances;
		vector<MATThresholds*> mats;
		double voltage, time;

};

#endif