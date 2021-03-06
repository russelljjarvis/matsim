# MATSim

MATSim is a package for simulation of the MAT (Multi-timescale Adaptive Threshold) model [(Kobayashi et al., 2009)][mat link], i.e., a non-resetting leaky membrane equipped with a dynamic firing threshold.

Two following modes of membrane stimulation are possible:
* Conductances (with reversal potentials) modelled as Poisson shot noise with an exponential envelope
* Conductances (with reversal potentials) modelled as an Ohrstein-Uhlenbeck process

![visualization of simulation](https://github.com/Tom83B/matsim/blob/master/examples/intro.png)

The package is written in C++ and wrapped with Cython for high efficiency.

## Installation

### Dependencies

MATSim requires:
* Python (>= 3.6)
* NumPy (>= 1.16.3)
* Pandas (>= 0.24.2)

Although older versions might be sufficient, they haven't been tested.

### User installation

#### Using pip

Install the package using `pip`:
```
pip install matsim
```

#### Compilation from source

First, use Cython to create the file `matsim.cpp` (in the directory `matsim`):
```
cython --cplus matsim.pyx
```
Then use the setup file to install the package:
```
python setup.py install
```

[mat link]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2722979/
