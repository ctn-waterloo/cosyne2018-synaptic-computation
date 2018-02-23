# Conductance-Based Synapses as a Computational Resource in the Neural Engineering Framework
## Andreas Stöckel, Aaron R. Voelker, Chris Eliasmith; presented at Cosyne 2018 in Denver, CO

![Picture with excerpts from the poster](/doc/header.png)

This repository accompanies our submission to Cosyne 2018 and demonstrates
synaptic computation in the NEF by exploiting nonlinear interactions of
conductance-based synapses.

## Abstract

Nonlinear interaction in the dendritic tree is known to be an important computational resource in biological neurons. Yet, high-level neural compilers ‒ such as the Neural Engineering Framework (NEF), or the predictive coding method published by Denève et al. in 2013 ‒ tend not to include conductance-based nonlinear synaptic interactions in their models, and so do not exploit these interactions systematically. In this study, we extend the NEF to include synaptic computation of nonlinear multivariate functions, such as controlled shunting, multiplication, and the Euclidean norm. We present a theoretical framework that provides sufficient conditions under which nonlinear synaptic interaction yields a similar precision compared to traditional NEF methods, while reducing the number of layers, neurons, and latency in the network. The proposed method lends itself to increasing the computational power of neuromorphic hardware systems and improves the NEF's biological plausibility by mitigating one of its long-standing limitations, namely its reliance on linear, current-based synapses. We perform a series of numerical experiments with a conductance-based two-compartment LIF neuron model. Preliminary results show that nonlinear interactions in conductance-based synapses are sufficient to compute a wide variety of nonlinear functions with performance competitive to using an additional layer of neurons as a nonlinearity.

* **Full Abstract** A PDF containing the full abstract submitted to Cosyne can be found [here](https://raw.githubusercontent.com/ctn-waterloo/cosyne2018-synaptic-computation/master/doc/astoeckel_cosyne_2018_abstract.pdf). *Note*: the preliminary results shown in the abstract are *obsolete* since we've switched to a QP solver instead of using gradient descent. The poster contains the most recent values for spiking simulations.
* **Poster** A low resolution [PDF copy](https://raw.githubusercontent.com/ctn-waterloo/cosyne2018-synaptic-computation/master/doc/astoeckel_cosyne_2018_poster_low_quality.pdf) of the poster presented at Cosyne 2018 can be found in the `doc` subdirectory.

## Running and Dependencies

Code is located in the Python 3 Jupyter notebook `synaptic_computation.ipynb`, as well as the `nef_synaptic_computation` package. A bare-bones version of the `nef_synaptic_computation` package is included with this repository. You can install it by running
```bash
pip3 install --user -e .
```
from the root directory of this repository (depending on your Python 3 distribution you may have to replace `pip3` with `pip`). In addition to Python 3.6 and a recent version of Jupyter you'll need to install recent versions of the following `pip` packages:

* numpy
* matplotlib
* nengo
* scipy
* cvxpy
* pandas

Note that spiking simulations will run significantly faster (about 100 times) on Linux machines if you have a recent version of the GNU C++ compiler installed.
