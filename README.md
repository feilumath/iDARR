# iDARR: Iterative Data-Adaptive RKHS Regularization

## Introduction

iDARR is an iterative regularization method tailored for large-scale inverse problems in the form `Ax = b`. It leverages the norm of a data-adaptive Reproducing Kernel Hilbert Space (RKHS) and is particularly effective when there is little prior information about the regularization norm for `x`.

The method searches for solutions in subspaces where the true solution can be identified, with the data-adaptive RKHS penalizing the spaces of small singular values. At the core of the method is a new generalized Golub-Kahan bidiagonalization procedure that recursively constructs orthonormal bases for a sequence of RKHS-restricted Krylov subspaces. The method is scalable with a complexity of $O(kmn)$ for $m$-by-$n$ matrices with $k$ denoting the iteration numbers.

For details of the method, refer to the paper: [Scalable Iterative Data-adaptive RKHS Regularization](https://arxiv.org/abs/2401.00656) by Haibo Li, Jinchao Feng, and Fei Lu.

## Features

iDARR is designed to:

- Provide an efficient and scalable solution to large-scale inverse problems.
- Adaptively determine the regularization norm based on data and the underlying model.
- Improve accuracy and stability in solutions of inverse problems with limited prior information.



## Usage

This repo implements iDARR in MATLAB. It uses the L-curve method for early stopping. Tests include the Frodholm integral equation and 2D image deblurring. The tests compares the DA-RKHS norm with two widely used L2 and l2 norms, denoted by IR-L2 and IR-l2 respectively. The L-curve method and the image handling tools in nice [IR Tools](https://link.springer.com/article/10.1007/s11075-018-0570-7) by Silvia Gazzola, Per Christian Hansen, and James G. Nagy.

- iDARR_demo.m:  This script provides a quick demonstration on using iDARR to solve discrete problems in the form `Ax = b`. 

- runme_main.m: This script contains the numerical tests on the Fredhold integral equation. It examines the stability of the method in different noise levels. The tests show that iDARR outperforms IR-L2 and IR-l2.    

- test_deblur.m: This script applies iDARR to 2D image deblurring, showcasing its potential in practical applications. 


---

For any issues or suggestions, feel free to contact us or open an issue on this repository.
