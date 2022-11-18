# Enumerating Grid Classes of Signed Permutations

An implementation of the algorithm described in [this paper](https://www.google.com/).  *link not set up yet since paper is pending publication*

This repository can be cloned into the working directory from the terminal using:

    git clone https://github.com/skora7/SignedPermutationClasses.git
    
### Prerequisites
This code works with Python 3.

In order to compute generating functions, either:
1. [SymPy](https://www.sympy.org/en/index.html) must be installed
2. The code must be run inside a [Sage](https://www.sagemath.org/) session
3. The FUNCTION_STRINGS constant on line 45 must be set to *True*

Everything else (including enumerating polynomials) will work fine with no additional dependencies.


### How To Use
The algorithm is implemented in permclass.py, which consists of three classes:
1. GridClass
2. SignedPerm
3. Polynomial

The purpose of each class is self-explanatory.  GridClass encapsulates all three major steps of our algorithm.  The **completion** and **compacting** steps are performed upon instantiation, while the **enumeration step** is performed any time *genfcn()* or *polynomial()* are called.
