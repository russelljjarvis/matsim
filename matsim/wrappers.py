from matsim import MATThresholds, Neuron

class LIF:
	def __init__(self, V_r=-80, R=50, C=0.1, V_thr=-60):
		self.V_r = V_r
		self.R = R
		self.C = C
		self.V_thr = V_thr

		self.threshold = MATThresholds(
		    alpha1=0,
		    alpha2=0,
		    tau1=10,
		    tau2=200,
		    omega=V_thr,
		    refractory_period=2,
		    resetting=True,
		    name='constant')

		self.neuron = Neuron(
		    resting_potential=V_r,
		    membrane_resistance=R,
		    membrane_capacitance=C,
		    mats=[self.threshold])

	def append_conductance(self, cond):
		self.neuron.append_conductance(cond)

	def timestep(self, dt):
		self.neuron.timestep(dt)

	@property
	def voltage(self):
		return self.neuron.voltage

	@property
	def time(self):
		return self.neuron.time

	def __getitem__(self, arg):
		if arg == 'V':
			return self.voltage

		elif arg == 't':
			return self.time

		else:
			pass