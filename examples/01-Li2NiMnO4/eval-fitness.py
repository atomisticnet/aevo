#!/usr/bin/env python

"""
Example fitness evaluation script for the aevo package.  In this
example, the electrostatic energy of an ionic structure is evaluated
using pymatgen.

For each file named 'POSCAR-*.vasp', this script generated an output
file with the electrostatic energy named 'POSCAR-*.vasp.energy'.

This script does not use command line options, so that the line where it
is called from the driver script does not have to be changed.  Instead,
modify all parameters directly within this script.

Adjustable parameters:
  oxidation_states (dict): dictionary with the oxidation states of all
    chemical species that can occur in the POSCAR files
  num_cores (int): number of cores for parallelization over files

"""

import argparse
import functools
import glob
import multiprocessing as mp

from pymatgen.io.vasp.inputs import Poscar
from pymatgen.analysis.ewald import EwaldSummation


__author__ = "Alexander Urban"
__email__ = "aurban@atomistic.net"
__date__ = "2020-05-29"
__version__ = "0.1"

oxidation_states = {"O": -2, "Li": 1, "Ni": 3.0, "Mn": 4}
num_cores = 1
poscar_files = glob.glob("POSCAR-*.vasp")


def compute_energy(structure_files, oxi, num_cores=1):
    eval_structure = functools.partial(_eval_structure, oxi)
    if num_cores == 1:
        results = [eval_structure(fname) for fname in structure_files]
    else:
        pool = mp.Pool(processes=num_cores)
        results = pool.map(eval_structure, structure_files)
    for filename, energy in results:
        idx = filename.replace("POSCAR-", "").replace(".vasp", "")
        with open("FITNESS.{}".format(idx), "w") as fp:
            fp.write("{}".format(energy))


def _eval_structure(oxi, filename):
    struc = Poscar.from_file(filename).structure
    struc.add_oxidation_state_by_element(oxi)
    assert(struc.charge == 0.0)
    energy = EwaldSummation(struc).total_energy
    return (filename, energy)


def get_fitness():
    compute_energy(poscar_files, oxi=oxidation_states, num_cores=num_cores)


if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(
        description=__doc__+"\n{} {}".format(__date__, __author__),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    args = parser.parse_args()
    get_fitness()
