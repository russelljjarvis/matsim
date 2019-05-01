# MATSim

MATSim is a package for simulation of the MAT (Multi-timescale Adaptive Threshold) model [(Kobayashi, 2009)][mat link], i.e., a non-resetting leaky membrane equipped with a dynamic firing threshold.

Two following modes of membrane stimulation are possible:
* Conductances (with reversal potentials) modelled as Poisson shot noise with an exponential envelope
* Conductances (with reversal potentials) modelled as an Ohrstein-Uhlenbeck process

The package is written in C++ and wrapped with Cython for high efficiency.

## Installation

Install the package using `pip`:
```
pip install matsim
```

[mat link]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2722979/