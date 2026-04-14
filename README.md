# MATLAB Implementation of DPIM

## Overview

This repository provides a code package for the **Direct Parametrisation of Invariant Manifolds (DPIM)**. It is mainly intended to demonstrate how DPIM can be implemented in practice with **MATLAB language** to compute:

- the **autonomous part** of **invariant manifold** of a nonlinear mechanical system
- and its **reduced-order model**.

The main purpose of this code is therefore **pedagogical and illustrative**: it helps users get started with MATLAB-based DPIM computation in a clear and accessible way.

## Scope

This MATLAB code focuses on the core ingredients required to obtain:

- a **high-order parametrisation of invariant manifolds**, and
- the corresponding **reduced dynamics** on the manifold.

It is particularly suitable for users who want to:

- learn the main workflow of DPIM,
- understand how invariant-manifold-based model reduction is implemented numerically,
- reproduce the essential reduced-order formulation in MATLAB.

## Limitations

This repository is **not intended to be a full-featured implementation** of all DPIM capabilities.

In particular, the code here mainly addresses the **autonomous part** of the reduced dynamics.  
For more advanced treatments, especially the **non-autonomous formulation**, users are encouraged to refer to the **Julia version of MORFE**, which provides a more complete and advanced implementation environment.

## References

The detailed theoretical background and computational principles can be found in the following references:

1. **Vizzaccaro, A., Opreni, A., Salles, L., et al.**  
   *High order direct parametrisation of invariant manifolds for model order reduction of finite element structures: application to large amplitude vibrations and uncovering of a folding point.*  
   **Nonlinear Dynamics**, 2022, **110**(1): 525-571.

2. **Wang, T., Touzé, C., Li, H. Q., et al.**  
   *Enhancing damping performance below the bandgap in metamaterial beams with geometric nonlinearity and bistable attachments via nonlinear energy transfer.*  
   **Mechanical Systems and Signal Processing**, 2025, **240**: 113362.

## Suggested Citation

If you use the ideas or implementation logic in this repository, please cite the references above as appropriate.
