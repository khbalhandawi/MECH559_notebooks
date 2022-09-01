### MECH 559 examples

## Prerequisites
- MATLAB (2019) or later
	- With the ``Optimization Toolbox``

- Julia 1.7 or later
	- IJulia
	- WebIO
	- Interact
	- LaTeXStrings
	- Plots
	- LatinHypercubeSampling
	- MLBase
	- Surrogates

### Installing Julia

The Julia interpreter and compiler can be downloaded from [http://julialang.org/downloads/](http://julialang.org/downloads/)  

Julia can be installed on major OS: Windows, MacOS X, Linux. Under Windows, it is possible to launch Julia via the start menu or, if present, by clicking on the appropriate icon. Under Linux, to launch the Julia interpreter in interactive mode, just enter at the terminal level

```
Julia
```

**Installing packages**

use the following command to install all the prerequisite packages for this course

```
Using Pkg
Pkg.add.(["IJulia", "WebIO", "Interact", "LaTeXStrings", "Plots", "LatinHypercubeSampling", "MLBase", "Surrogates"])
```

## List of examples
- [Multi-objective optimization of beam example](./Beam_example/) 
- [Data fits using surrogate models](./data_fits/) 
- Unconstrained optimization ([gradient descent](./Unc/) , [Rosenbrock](./Rosenbrock/))
- [Linear programming example](./LP/) 
- [Reduced gradient visualization](<./Reduced Gradient/>) 
- [Augmented Lagrangian Example](./AUgmentedLagrangianExample/) 
- [Nonlinear programming examples](./NLPexample/) 