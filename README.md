===================================
# PYTHON-ISOTOPIC-DEPLETION-PACKAGE
===================================

Custom depletion tool for transmutation and decay heat analysis.

Description
============
The solver can be used to solve the following problems:
1. Depletion (transmutation + decay)
2. Only decay
3. Variable power or flux history
4. Allows to define cross section sets for different state points.

The package calculates the following time-dependent results:
- Nuclide densities,
- Activities,
- Decay heat,
- Radiotoxicity,

The Bateman solutions are solved using the following methods:
1. The Chebyshev Rational Approximation Method (CRAM) 
2. The expm built-in function that relies on the Pade approximation.
3. Built-in odeint solver

Installation
============
After downloading all the package navigate to the root directory.
The ``setup.py`` defines how to install the code.
Execute the following line:

NO NEED:: python setup.py sdist bdist_wheel

python setup.py install --user

The command should create a ``dist`` directory where you should find the `<tar.gz>`
and `<.whl>` files. The `<tar.gz>` file is a source archive and `<.whl>` is a built distribution.

Before installing you may need to make sure that you have the latest pip, setuptools,and wheel on your system.
Follow the next execution line:

`<python -m pip install -upgrade pip setuptools wheel>`


Execution
=========
The `pyIsoDep` package can now be imported as as standard PyPI using:

import `import pyIsoDep`

Currently a memory-based input data is required.
The user must import the following modules:

`from pyIsoDep.functions.generatedata TransmutationData


References
==========
[1] Higham, N. J., “The Scaling and Squaring Method for the Matrix Exponential Revisited,” SIAM J. Matrix Anal. Appl., 26(4) (2005), pp. 1179–1193.

[2] Al-Mohy, A. H. and N. J. Higham, “A new scaling and squaring algorithm for the matrix exponential,” SIAM J. Matrix Anal. Appl., 31(3) (2009), pp. 970–989.

[3] Golub, G. H. and C. F. Van Loan, Matrix Computation, p. 384, Johns Hopkins University Press, 1983.

[4] Moler, C. B. and C. F. Van Loan, “Nineteen Dubious Ways to Compute the Exponential of a Matrix,” SIAM Review 20, 1978, pp. 801–836. Reprinted and updated as “Nineteen Dubious Ways to Compute the Exponential of a Matrix, Twenty-Five Years Later,” SIAM Review 45, 2003, pp. 3–49.
