# Cosyne 2018 Project: Synaptic Computation in the NEF
## Andreas Stöckel, Aaron R. Voelker, Chris Eliasmith

This repository accompanies our submission to Cosyne 2018 and demonstrates
synaptic computation in the NEF by exploiting nonlinear interactions of
conductance-based synapses.

## Abstract

Nonlinear interaction in the dendritic tree is known to be an important computational resource in biological neurons. Yet, high-level neural compilers -- such as the Neural Engineering Framework (NEF), or the predictive coding method published by Denève et~al. in 2013 -- tend not to include conductance-based nonlinear synaptic interactions in their models, and so do not exploit these interactions systematically. In this study, we extend the NEF to include synaptic computation of nonlinear multivariate functions, such as controlled shunting, multiplication, and the Euclidean norm. We present a theoretical framework that provides sufficient conditions under which nonlinear synaptic interaction yields a similar precision compared to traditional NEF methods, while reducing the number of layers, neurons, and latency in the network. The proposed method lends itself to increasing the computational power of neuromorphic hardware systems and improves the NEF's biological plausibility by mitigating one of its long-standing limitations, namely its reliance on linear, current-based synapses. We perform a series of numerical experiments with a conductance-based two-compartment LIF neuron model. Preliminary results show that nonlinear interactions in conductance-based synapses are sufficient to compute a wide variety of nonlinear functions with performance competitive to using an additional layer of neurons as a nonlinearity.

## Running and Dependencies

All code is located in the Python 3 Jupyter notebook `synaptic_computation.ipynb`. In addition to Python 3.6 and a recent version of Jupyter you'll need to install recent versions of the following `pip` packages:

* numpy
* matplotlib
* nengo
* scipy
* pandas
