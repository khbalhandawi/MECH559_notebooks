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

### Installing Jupyter

To install Jupyter you need to have an existing Python installation (Anaconda or otherwise). You may then use the following commands in the console

* Linux/MacOS
```
pip3 install jupyter
python3 -m pip install --upgrade webio_jupyter_extension
```

* Windows
```
pip install jupyter
python -m pip install --upgrade webio_jupyter_extension
```

### Installing Julia

The Julia interpreter and compiler can be downloaded from [http://julialang.org/downloads/](http://julialang.org/downloads/)  

Julia can be installed on major OS: Windows, MacOS X, Linux. Under Windows, it is possible to launch Julia via the start menu or, if present, by clicking on the appropriate icon. Under Linux, to launch the Julia interpreter in interactive mode, just enter at the terminal level

```
Julia
```

Note for windows: Be sure to check the "add to PATH" option during installation to add Julia to your environment variables.

**Installing packages**

use the following command to install all the prerequisite packages for this course

```
using Pkg
Pkg.add.(["IJulia", "WebIO", "Interact", "LaTeXStrings", "Plots", "LatinHypercubeSampling", "MLBase", "Surrogates"])
```

## List of examples
- [Multi-objective optimization of beam example](./1_Beam_example/) 
- [Monotonicity theorem examples](./2_monotonicity/) 
- [Data fits using surrogate models](./3_data_fits/) 
- Unconstrained optimization ([gradient descent](./Unc/) , [Rosenbrock](./Rosenbrock/))
- [Linear programming example](./LP/) 
- [Reduced gradient visualization](<./Reduced Gradient/>) 
- [Augmented Lagrangian Example](./AUgmentedLagrangianExample/) 
- [Nonlinear programming examples](./NLPexample/) 