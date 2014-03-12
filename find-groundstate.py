#!/usr/bin/env python

"""
Search for the ground state ordering of a crystal structure with
multiple sub-lattices, using an evolutionary algorithm.
"""

__author__ = "Alexander Urban"
__date__   = "2014-03-03"

import shutil
import argparse

import pygs

#----------------------------------------------------------------------#

def find_groundstate(infile, statefile=None, outdir=None, print_top_N=None):

    if statefile is None:
        print
        print " Generating new population of trials."
        print
        ga = pygs.Evolution.from_parameter_file(infile)
        ga.write_unevaluated(directory=outdir)
        statefile = 'state.json'
        print " Saving state to: {}".format(statefile)
        with open(statefile, 'w') as fp:
            fp.write(ga.to_JSON())
        print
    else: # There is a state file.  Read it.
        print
        print " Reading state file: {}".format(statefile)
        with open(statefile, 'r') as fp:
            ga = pygs.Evolution.from_JSON(fp)
        ga.update_parameters(infile)
        for trial in ga.trials: # (fixme: remove this loop later)
            if (trial.id not in ga.fitness_history
               ) and (trial.fitness is not None):
                ga.fitness_history[trial.id] = trial.fitness
        if print_top_N is not None:
            # just print the top N trials, no GA
            print
            print " Top {} trials after {} generations:".format(
                print_top_N, ga.generation-1)
            print
            for (ID,fitness) in ga.top_N(print_top_N):
                print "   {}  {:15.8e}  !".format(ID, fitness)
            print
        else:
            # actual GA
            shutil.move(statefile, statefile + '-bak')
            ga.read_fitness(directory=outdir)
            ga.select()
            print (" Generation   Average fitness   Minimum fitness   "
                  + "Maximum fitness")
            print " {:10d}   {:15.8e}   {:15.8e}   {:15.8e}   !".format(
                ga.generation, ga.average_fitness, ga.min_fitness,
                ga.max_fitness)
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
        help    = "Path to a state file containing restart information "
                  + "(JSON format).",
        default = None)

    parser.add_argument(
        "--directory", "-d",
        help    = "Path to directory for input/output files.  Per default the "
                  + "directory name is the current generation.",
        default = None)

    parser.add_argument(
        "--top-N",
        help    = "Print the IDs of the best N trials and exit.",
        type    = int,
        default = None)

    args = parser.parse_args()

    find_groundstate(args.input_file,
                     statefile=args.state,
                     outdir=args.directory,
                     print_top_N=args.top_N)
