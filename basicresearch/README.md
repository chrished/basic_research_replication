Back to Basics
==============

Replication code for
```
Back to Basics: Basic Research Spillovers, Innovation Policy and Growth by Ufuk Akcigit, Douglas Hanley, and Nicolas Serrano-Velarde. Review of Economic Studies, 2020.
```

Many configuration options are set in `inialg.m`. This includes which parameter values to use (from the `params` directory) and which equilibrium values to use as a starting guess (from the `eqvars` folder) as well as moment targets in `targets`. In addition, there are a number of algorithmic parameters and certain special policy related flags set there.

The entry point for solving the model is `eqsolve.m`. Note that you must run `eqsolve(1)` to save the values required for future consumption equivalent calculations. The main estimation can be run from `smm.m`. Here you can choose which parameters to fix and which to optimize over by index. Various types of optimal policy can be computed from `optpol.m`, and the social planner's problem is solved using `socplan.m`.

The system of equations characterizing the model solution is computed by `eqfunc.m`. Drilling down through the various routines called from this file will give you a good idea of the approach used. Much of the actual computational time is spent in `qdistbin.m` finding the distributions over `qhat` using a null space approach. The other major offender is `nmdist.m` which computes the joint distribution over the number of industries `m` and the number of products `n` using an iterative approach. To compute the necessary moments, one actually needs the joint distribution over `n`, `m`, and the portfolio of `qhat` values. This requires a more detailed firm-level simulation which is implemented in `compfs.m`.

We do make use of the `mex` subsystem to incorporate C++ code. These files can be found in the `standalone` directory. First, the iterative scheme for finding the joint `m`/`n` distribution can be found in `eigsim`. Second, the detailed firm simulation can be found in `firmsim`. Both of these have `cpp` and `cu` implementations corresponding to the CPU and GPU versions, respectively. We recommend using the CPU implementation as the speedup from using the GPU ended up being minimal and in some cases negative. To compile the CPU version, simply run `make` in the `standalone` directory and symbolically link or copy the resulting `.mexa64` files to the parent directory.
