#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <ctime>
#include "bio.h"
using namespace std;

#define PI           3.14159265358979323846  /* pi */

int poisson_rand(double rate) {
	double randu = ((double) rand() / (RAND_MAX));
	double cummulation = 0;
	double pmf = exp(-rate);
	int poisson_rand = 0;

	// cout << endl << "randu: " << randu << endl;

	while (randu > cummulation) {
		// cout << "cummulation" << cummulation << endl;
		cummulation += pmf;
		poisson_rand += 1;
	    pmf *= rate / poisson_rand;
	}
	// cout << "finished with   " << cummulation << endl;

	return poisson_rand - 1;
}

double normal_rand() {
	double randu1 = ((double) rand() / (RAND_MAX));
	double randu2 = ((double) rand() / (RAND_MAX));

	// cout << randu1 << endl << randu2 << endl << log(randu1) << endl;
	// cout << cos(2 * PI * randu2) << endl;
	// cout << -2 * log(randu1) * cos(2 * PI * randu2) << endl << endl;

	return sqrt(-2 * log(randu1)) * cos(2 * PI * randu2);
}

void Conductance::update(double dt) { }
void Conductance::set_rate(double rate) { }
double Conductance::get_g() {return 0;}
double Conductance::get_reversal() {return 0;}
void Conductance::activate() { }

ExponentialConductance::ExponentialConductance() {}
double ExponentialConductance::get_g() {return this->g;}
double ExponentialConductance::get_reversal() {return this->reversal;}

ExponentialConductance::ExponentialConductance(double g_peak, double reversal, double decay) {
	this->g_peak = g_peak;
	this->reversal = reversal;
	this->decay = decay;
	this->g = 0;
}

void ExponentialConductance::update(double dt) {
	g *= exp(-dt / decay);
}

void ExponentialConductance::activate() {
	g += g_peak;
}

ShotNoiseConductance::ShotNoiseConductance() {}
double ShotNoiseConductance::get_g() {return this->g;}
double ShotNoiseConductance::get_reversal() {return this->reversal;}
void ShotNoiseConductance::activate() { }

ShotNoiseConductance::ShotNoiseConductance(double rate, double g_peak, double reversal, double decay) {
	srand(std::time(nullptr));

	this->rate = rate;
	this->g_peak = g_peak;
	this->reversal = reversal;
	this->decay = decay;
	this->g = 0;
}

void ShotNoiseConductance::update(double dt) {
	int n_spikes = poisson_rand(rate * dt);
	g *= exp(-dt / decay);

	g += n_spikes * g_peak * decay * (1 - exp(-dt / decay)) / dt;
}

void ShotNoiseConductance::set_rate(double rate) {
	this->rate = rate;
}

OUConductance::OUConductance() {}
double OUConductance::get_g() {return this->g;}
double OUConductance::get_reversal() {return this->reversal;}
void OUConductance::activate() { }

OUConductance::OUConductance(double rate, double g_peak, double reversal, double decay) {
	srand(std::time(nullptr));

	this->rate = rate;
	this->g_peak = g_peak;
	this->reversal = reversal;
	this->decay = decay;

	this->set_rate(rate);
	this->g = this->mean;
	// cout << normal_rand() << endl;
	// cout << g << endl;
}

double OUConductance::get_A(double dt) {
	double A = sqrt(D * decay / 2 * (1 - exp(-2 * dt / decay)));
	return A;
}

void OUConductance::update(double dt) {
    double A = this->get_A(dt);

    // cout << g << "   " << mean << endl;
    g = mean + (g - mean) * exp(-dt / decay) + A * normal_rand();
}

void OUConductance::set_rate(double rate) {
	// cout << "funguju?";
	this->mean = decay * g_peak * rate;
	this->sigma = sqrt(decay * rate / 2) * g_peak;
	this->D = 2 * (sigma * sigma) / decay;
}

MATThresholds::MATThresholds(double alpha1, double alpha2, double tau1, double tau2,\
	double omega, double refractory_period, bool resetting = false) {
	this->alpha1 = alpha1;
	this->alpha2 = alpha2;
	this->tau1 = tau1;
	this->tau2 = tau2;
	this->omega = omega;
	this->refractory_period = refractory_period;
	this->resetting = resetting;

	t1 = 0;
	t2 = 0;
	threshold = t1 + t2 + omega;
	past_spike_time = refractory_period;
}

void MATThresholds::fire(double time) {
	if (past_spike_time >= refractory_period) {
		t1 += alpha1;
		t2 += alpha2;
		spike_times.push_back(time);
		past_spike_time = 0;
	}
}

void MATThresholds::update(double dt) {
	past_spike_time += dt;
	t1 *= exp(-dt / tau1);
	t2 *= exp(-dt / tau2);
	threshold = t1 + t2 + omega;
}

vector<double> MATThresholds::get_spike_times() {
	return spike_times;
}

void MATThresholds::reset_spike_times() {
	vector<double> spike_times;
	this->spike_times = spike_times;
}

Neuron::Neuron() {}

Neuron::Neuron(double resting_potential, double membrane_resistance, double membrane_capacitance, vector<MATThresholds*> mats) {
	this->resting_potential = resting_potential;
	this->membrane_resistance = membrane_resistance;
	this->membrane_capacitance = membrane_capacitance;
	this->mats = mats;

	voltage = resting_potential;
	time_constant = membrane_capacitance * membrane_resistance;
	time = 0;
}

void Neuron::append_conductance(Conductance* conductance) {
	conductances.push_back(conductance);
}

void Neuron::integrate_voltage(double dt) {
	double tot_conductance = 0;
	double tot_gr = 0;
	double factor, v_bar, tau;

	for (auto c : conductances) {
		tot_conductance += c->get_g();
		tot_gr += c->get_g() * c->get_reversal();
	}

	v_bar = (resting_potential + membrane_resistance * tot_gr) / (1 + membrane_resistance * tot_conductance);
	tau = time_constant / (1 + membrane_resistance * tot_conductance);

	voltage = v_bar + (voltage - v_bar) * exp(-dt / tau);
}

void Neuron::timestep(double dt) {
	time += dt;

	for (auto c : conductances) {
		c->update(dt);
	}

	this->integrate_voltage(dt);

	for (auto mat : mats) {
		mat->update(dt);
		if (mat->threshold <= voltage) {
			mat->fire(time);
			
			if (mat->resetting == true) {
				this->voltage = this->resting_potential;
			}

			for (auto c : conductances) {
				c->activate();
			}
		}
	}
}

HHNeuron::HHNeuron() {
	double A = 1e-4;

	this->g_l = 0.045 * 1000 * A;
	this->E_l = -80;
	this->c_m = 1 * A * 1000.;
	this->E_na = 50;
	this->g_na = 50 * 1000. * A;
	this->E_k = -90;
	this->g_k = 5 * 1000 * A;

	this->V = E_l;
	this->time = 0;

	this->m = 0;
	this->h = 0;
	this->n = 0;

	this->i_na = 0;
	this->i_k = 0;
	this->i_l = 0;
}

void HHNeuron::append_conductance(Conductance* conductance) {
	conductances.push_back(conductance);
}

void HHNeuron::integrate_voltage(double dt) {
	double i_syn = 0;
	double am, ah, an, bm, bh, bn;
	double dVdt;
	
	double VT = -58;
	double VS = -10;

	for (auto c : conductances) {
		i_syn += c->get_g() * (c->get_reversal() - V);
	}
	
	am = -0.32 * (V - VT - 13) / (exp(-(V - VT - 13) / 4) - 1);
	bm = 0.28 * (V - VT - 40) / (exp((V - VT - 40) / 5) - 1);

	ah = 0.128 * exp(-(V - VT - VS - 17) / 18);
	bh = 4 / (1 + exp(-(V - VT - VS - 40) / 5));

	an = -0.032 * (V - VT - 15) / (exp(-(V - VT - 15) / 5) - 1);
	bn = 0.5 * exp(-(V - VT - 10) / 40);
	// am = (0.182 * (V + 35)) / (1 - exp(-(V + 35) / 9));
	// bm = -(0.124 * (V + 35)) / (1 - exp((V + 35) / 9));

	// ah = 0.25 * exp(-(V + 90.) / 12.);
	// bh = 0.25 * exp((V + 62.) / 6.) / exp((V + 90.) / 12.);

	// an = (0.02 * (V - 25.) / 9.) / (1 - exp(-(V - 25.) / 9.));
	// bn = -0.002 * (V - 25.) / (1 - exp((V - 25.) / 9));

	this->i_na = g_na * m * m * m * h * (V - E_na);
	this->i_k  = g_k  * n * n * n * n * (V - E_k);
	this->i_l  = g_l * (V - E_l);

	dVdt = (-this->i_l - this->i_na - this->i_k + i_syn) / c_m;
	this->V += dVdt * dt;
	this->m += dt * (am * (1 - m) - bm * m);
	this->h += dt * (ah * (1 - h) - bh * h);
	this->n += dt * (an * (1 - n) - bn * n);

	// cout << endl;
	// cout << "i_na:   " << i_na << endl;
	// cout << "i_k:   " << i_k << endl;
	// cout << "i_l:   " << i_l << endl;
	// cout << "i_s:   " << i_syn << endl;
	// cout << "V:   " << V << endl;
	// cout << endl;

	// cout << endl;
	// cout << "m:   " << m << endl;
	// cout << "h:   " << h << endl;
	// cout << "n:   " << n << endl;
	// cout << endl;
}

void HHNeuron::timestep(double dt) {
	time += dt;

	for (auto c : conductances) {
		c->update(dt);
	}

	this->integrate_voltage(dt);
}

// int main(int argc, char const *argv[])
// {
// 	srand(time(nullptr));

// 	HHNeuron neuron;
// 	ShotNoiseConductance exc (10, 0.0015, 0, 3);
// 	ShotNoiseConductance inh (5, 0.0015, -75, 10);

// 	neuron.append_conductance(&inh);
// 	neuron.append_conductance(&exc);

// 	for (int i=0; i<10000; i++) {
// 		neuron.timestep(0.025);
// 		cout << i << endl;
// 		// if (i % 10000 == 0) cout << i << endl;
// 		// cout << neuron.voltage << "   " << neuron.mats[0]->threshold << endl;
// 	}

// 	// cout << mat.get_spike_times().size();

// 	return 0;
// }
