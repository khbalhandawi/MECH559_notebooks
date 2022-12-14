# MECH 559 examples
## List of examples
- [Multi-objective optimization of beam example](./1_Beam_example/) 
- [Monotonicity and extreme value theorem examples](./2_monotonicity_boundedness/) 
- [Data fits using surrogate models](./3_data_fits/) 
- [Concepts and examples for unconstrained optimization (FONCs, SOSCs)](./4_unconstrained/) 
- [Unconstrained optimization algorithms (Gradient descent and variants)](./5_unconstrained_algorithms/)
- [Linear programming examples](./6_linear_programming/) 
- [Non linear programming examples (KKT conditions and regular points)](<./7_nonlinear_programming/>) 
- [Non linear programming algorithms and examples](./8_nonlinear_programming_algorithms/) 
- [Derivative free optimization examples](./9_dfo/) 
- [Multi-disciplinary optimization (MDO) and Multi-disciplinary analysis (MDA) examples](./10_mdo/) 
- [Source files for compiling a blackbox](./Blackbox_files/)

## Prerequisites
- MATLAB (2019) or later with
	- ``Optimization Toolbox``
	- ``Global Optimization Toolbox``
	- ``Symbolic Math Toolbox``

- Julia 1.7 or later with
	- IJulia
	- WebIO
	- Interact
	- LaTeXStrings
	- Plots
	- PyPlot
	- LatinHypercubeSampling
	- MLBase
	- Surrogates
	- CSV

- Python 3.8 or later with dependencies in [requirements.txt](./requirements.txt)

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
Pkg.add.(["IJulia", "WebIO", "Interact", "LaTeXStrings", "Plots", "PyPlot", "LatinHypercubeSampling", "MLBase", "Surrogates", "CSV"])
```

### Setting up your Python environment

**MacOS/Linux**

Create a virtual environment where you can install the required packages for this class locally

```
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

**Windows**
```
python -m venv .env
.env\Scripts\activate
pip install -r requirements.txt
```