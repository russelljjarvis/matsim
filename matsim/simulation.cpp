#include <vector>
#include <iostream>
#include <cmath>
#include "simulation.h"
#include "bio.h"

using namespace std;

Simulation::Simulation(Neuron neuron, ShotNoiseConductance exc, ShotNoiseConductance inh) {
	this->neuron = neuron;
	this->exc = exc;
	this->inh = inh;

	this->neuron.append_conductance(&(this->exc));
	this->neuron.append_conductance(&(this->inh));
}

void Simulation::run(double time, double dt, double exc_rate, double inh_rate) {
	(*neuron.conductances[0]).set_rate(exc_rate);
	(*neuron.conductances[1]).set_rate(inh_rate);

	double sim_time = 0;

	while (sim_time <= time) {
		neuron.timestep(dt);
		sim_time += dt;
	}
}

vector<int> sr_experiment(Neuron neuron, double time_window, double dt,
		vector<double> exc_intensities, vector<double> inh_intensities, int seed) {
	srand(seed);
	vector<int> spike_counts;

	int len = exc_intensities.size();
	int print_step = len / 100;
	if (print_step == 0) print_step = 1;

	for (int i=0; i<len; i++) {
		if (i % print_step == 0) {
			cout << i << " / " << len << endl;
		}
		(*neuron.conductances[0]).set_rate(exc_intensities[i]);
		(*neuron.conductances[1]).set_rate(inh_intensities[i]);

		double sim_time = 0;

		while (sim_time <= time_window) {
			neuron.timestep(dt);
			sim_time += dt;
		}

		for (auto mat_ptr : neuron.mats) {
			spike_counts.push_back((mat_ptr->get_spike_times()).size());
			mat_ptr->reset_spike_times();
		}
	}

	return spike_counts;
}

vector<vector<double> > sr_experiment_spike_times(Neuron neuron, double time_window, double dt,
		vector<double> exc_intensities, vector<double> inh_intensities, int seed) {
	srand(seed);
	vector<vector<double> > spike_times;

	int len = exc_intensities.size();
	int print_step = len / 100;
	if (print_step == 0) print_step = 1;

	for (int i=0; i<len; i++) {
		neuron.time = 0;

		// if (i % print_step == 0) {
		// 	cout << i << " / " << len << endl;
		// }
		(*neuron.conductances[0]).set_rate(exc_intensities[i]);
		(*neuron.conductances[1]).set_rate(inh_intensities[i]);

		double sim_time = 0;

		while (sim_time <= time_window) {
			neuron.timestep(dt);
			sim_time += dt;
		}

		for (auto mat_ptr : neuron.mats) {
			spike_times.push_back(mat_ptr->get_spike_times());
			mat_ptr->reset_spike_times();
		}
	}

	return spike_times;
}

int main(int argc, char const *argv[])
{
	srand(time(nullptr));

	vector<int> spike_counts;
	vector<double> intensities;

	MATThresholds mat(10, 1, 10, 200, -60, 2);
	vector<MATThresholds*> mats;
	mats.push_back(&mat);

	Neuron neuron(-80, 50, 0.1, mats);
	ShotNoiseConductance exc (10, 0.0015, 0, 3);
	ShotNoiseConductance inh (10, 0.0015, -75, 10);

	neuron.append_conductance(&exc);
	neuron.append_conductance(&inh);

	vector<double> exc_intensities;
	vector<double> inh_intensities;

	for (int i=0; i<100000; i++) {
		int exc_intensity = 5;
		int inh_intensity = 7;

		exc_intensities.push_back(exc_intensity);
		inh_intensities.push_back(inh_intensity);
	}

	spike_counts = sr_experiment(neuron, 1000, 0.1, exc_intensities, inh_intensities, 1);

	// Simulation sim(neuron, exc, inh);

	// for (int i=0; i<100; i++) {
	// 	sim.run(1000, 0.1, 10, 10);

	// 	for (auto mat_ptr : sim.neuron.mats) {
	// 		spike_counts.push_back((mat_ptr->get_spike_times()).size());
	// 		mat_ptr->reset_spike_times();
	// 	}
	// }
	
	for (auto cnt : spike_counts) {
		cout << cnt << endl;
	}

	cout << "???";

	return 0;
}