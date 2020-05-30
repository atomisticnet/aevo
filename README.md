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

[1] *ævo* reference:  N. Artrith, A. Urban, and G. Ceder, [*J. Chem. Phys.* **148** (2018) 241711](https://doi.org/10.1063/1.5017661).<br/>
[2] *Differential evolution* method: R. Storn and K. Price, [*J. Global Optim.* **11** (1997) 341-359](https://doi.org/10.1023/A:1008202821328).<br/>
[3] *pymatgen* toolkit:  S. P. Ong et al., [*Comput. Mater. Sci.* **68** (2013) 314–319](https://doi.org/10.1016/j.commatsci.2012.10.028).

## Installation

Make sure that all requirements have been installed.  Then install *ævo*
with [`pip`](https://pip.pypa.io) from within the git repository with

     $ pip install . --user

Once installed, make sure that the command line tool `aevolution.py` is
available in the system path with

     $ aevolution.py --help

which should print out the following help text:

     usage: aevolution.py [-h] [--state STATE] [--directory DIRECTORY]
                          [--print-top-N PRINT_TOP_N] [--save-top-N SAVE_TOP_N]
                          [input_file]

     Search for the ground state ordering of a crystal structure with
     multiple sub-lattices, using an evolutionary algorithm.

     2014-03-03 Alexander Urban, Nongnuch Artrith

     positional arguments:
       input_file            Input file in JSON format.

     optional arguments:
       -h, --help            show this help message and exit
       --state STATE, -s STATE
                             Path to a state file containing restart information
                             (JSON format).
       --directory DIRECTORY, -d DIRECTORY
                             Path to directory for input/output files. Per default
                             the directory name is the current generation.
       --print-top-N PRINT_TOP_N
                             Print the IDs of the best N trials and exit.
       --save-top-N SAVE_TOP_N
                             Save the best N trials into a separate directory and
                             exit.

## Usage

Evolutionary optimization using *ævo* proceeds via the following steps

0. Generate initial trial atomic configurations
1. Evaluate the fitness of all new trials in current population using a
   user-provided script
2. Cross trials from current population to generate new trials,
   select current trials with good fitness with greater probability
3. Introduce random mutations in the newly generated trials
4. Select which trials will be kept in the population with a probability
   proportional to the fitness
5. Continue with step 1.

The *crossing* step 2 follows the *differential evolution* method [2].

### Principal input file

Adjustable parameters are defined in an input file in [JSON
format](https://www.json.org).  Here an example:

     {
         "size"   : 15,
         "keep"   : 4,
         "mutate" : 0.1,
         "sites"  : "LiNiO2-4x4x1.vasp",
         "atom_types" : [ "Li", "Ni", "Mn", "O" ],
         "sublattices" : {
             "A" : { "Li": 8 },
             "B" : { "Ni": 8, "Mn": 8 },
             "C" : { "O": 32 }
         },
         "initial_population" : []
     }

- The `size` keyword defines the population size of a new generation
  (number of new trials per iteration).
- `keep` defines how many of the best trials will be guaranteed to be
  kept in the population.  All other trials will be selected
  stochastically with a probability that is proportional to their
  fitness.
- `mutate` defines the *mutation rate*, i.e., the probability that a
  trial is randomly modified (1.0 = 100%).
- `sites` specifies the path to an atomic structure file in [VASPS
  POSCAR format](https://www.vasp.at/wiki/index.php/Input#POSCAR) that
  defines the coordinates of the lattice sites and assigns them to
  sublattices.
- `atom_types` is a list of all allowed atomic species.  Vacancies do
  not have to be specified as a separate species.
- `sublattices` is a dictionary that defines the composition of each
  sublattice.  If fewer atoms are specified than sites contained in the
  sublattice, the remaining sites are taken to be vacancies.  The labels
  of the sublattices (`A`, `B`, and `C` in the example) are defined in
  the `sites` file.
- `initial_population` is a list of POSCAR files to be used as initial
  population.  If no paths are given or if this parameter is absent,
  random atomic configurations will be generated as initial trial
  population.

### Site definitions

The `sites` file is a regular POSCAR file that defines sublattice labels
instead of chemical species.

Example header of a `sites` file:

     A1 B1 C2
     1.0
         11.1551200000000      0.0000000000000     -0.0000000000000
          5.5770480000000     -1.4232520000000    -10.2143000000000
          1.3937190000000      4.6309280000000     -1.2767880000000
     A  B  C
     16 16 32
     direct
     ...

Each sublattice can be occupied by atoms of multiple chemical species,
as defined in the principal input file described above.

### Driver script

A *driver script* is needed to connect *ævo* with the user-defined
fitness evaluation script.  Usage examples of *ævo* optimization runs
that also include driver scripts and fitness evaluation scripts are
collected in the [examples subdirectory](./examples/).
