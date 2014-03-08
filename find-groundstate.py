#!/usr/bin/env python

"""
Insert description here.
"""

__author__ = "Alexander Urban"
__date__   = "2014-03-03"

import shutil
import argparse

import pygs

#----------------------------------------------------------------------#

def find_groundstate(infile, statefile=None, outdir=None):

    if statefile is None:
        print
        print " Generating new population of trials."
        print
        ga = pygs.Evolution.from_parameter_file(infile)
        ga.write_unevaluated(directory=outdir)
        statefile = 'state.json'
    else:
        print
        print " Restarting from state file: {}".format(statefile)
        with open(statefile, 'r') as fp:
            ga = pygs.Evolution.from_JSON(fp)
        shutil.move(statefile, statefile + '-bak')
        ga.update_parameters(infile)
        ga.read_fitness(directory=outdir)
        ga.select()
        print " Generation   Average fitness   Minimum fitness   Maximum fitness"
        print " {:10d}   {:15.8e}   {:15.8e}   {:15.8e}   !".format(
            ga.generation, ga.average_fitness, ga.min_fitness, ga.max_fitness)
        print
        ga.mate()
        ga.write_unevaluated(directory=outdir)
        ga.write_current_best('current_best.vasp')

    print " Saving state to: {}".format(statefile)
    with open(statefile, 'w') as fp:
        fp.write(ga.to_JSON())

    print

#----------------------------------------------------------------------#

if (__name__ == "__main__"):

    parser = argparse.ArgumentParser(
        description     = __doc__+"\n{} {}".format(__date__,__author__),
        formatter_class = argparse.RawDescriptionHelpFormatter )

    parser.add_argument(
        "input_file",
        help    = "Input file in JSON format.")

    parser.add_argument(
        "--state", "-s",
        help    = "Path to a state file containing restart information.",
        default = None)

    parser.add_argument(
        "--directory", "-d",
        help    = "Path to directory for output files.",
        default = None)

    args = parser.parse_args()

    find_groundstate(args.input_file,
                     statefile=args.state,
                     outdir=args.directory)
