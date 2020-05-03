Unbiased Estimation of the Gradient of the Log-Likelihood in Inverse Problems
==
This is the offical implementation for our [paper](https://arxiv.org/abs/2003.04896).

Introduction
--
We consider the problem of estimating a parameter $\theta \in \Theta$ associated to a Bayesian inverse problem. Treating the unknown initial condition as a nuisance parameter, typically, one must resort to a numerical approximation of gradient of the log-likelihood and, in particular, adopt a discretization of the problem in space or time. We develop a new methodology to unbiasedly estimate the gradient of the log-likelihood w.r.t. the unknown parameter, that is, an estimate which loses the discretization bias in expectation. Such a property is not only useful w.r.t. estimation in terms of the original stochastic model of interest, but can be used in stochastic gradient algorithms which benefit from unbiased estimates. We prove, under assumptions, that our estimator is not only unbiased but of finite variance. In addition, in comparison to multilevel estimation methods, when implemented on a single processor, we show that the cost to achieve a given level of error is comparable, not only practically, but theoretically. However, the new algorithm is highly amenable to parallel computation.

Installation
--
Install the package smccbase14_0.1.0.tar.gz in R before running code.
