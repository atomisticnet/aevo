ævo – Atomistic Evolution
--------------------------------------

A Python package for the optimization of atomic configurations using an
evolutionary algorithm.

*ævo* implements the *differential evolution* method [1] to search for the
global ground state atomic ordering with respect to a fitness function.
The package is not limited to a specific fitness function but instead
calls an external user-provide script to evaluate the fitness of a given
atomic configuration.

## Prerequisites

The *ævo* package makes use of the excellend *Python Materials Genomics*
(*pymatgen*) toolkit for the input and output of atomic structures.  See
the *pymatgen* website (http://pymatgen.org) for installation
instructions.

## References

Please cite the following references in all products/publications that
made use of ævo:

*ævo* reference:  N. Artrith, A. Urban, G. Ceder, [*J. Chem. Phys.* **148** (2018) 241711](https://doi.org/10.1063/1.5017661).

*pymatgen* reference:  S. P. Ong et al., [*Comput. Mater. Sci.* **68** (2013) 314–319](doi:10.1016/j.commatsci.2012.10.028).

## Installation

Make sure that all requirements have been installed.  Then install *ævo*
with [~pip~](https://pip.pypa.io) from within the git repository with

     $ pip install . --user
