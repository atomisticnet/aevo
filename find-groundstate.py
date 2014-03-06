#!/usr/bin/env python

"""
Insert description here.
"""

__author__ = "Alexander Urban"
__date__   = "2014-03-03"

import argparse

import pygs

#----------------------------------------------------------------------#

def find_groundstate(infile, statefile=None):

    if statefile is None:
        print
        print " Generating new population of trials"
        ga = pygs.Evolution.from_parameter_file(infile)
        ga.write_unevaluated()
    else:
        print
        print " Restarting from state file {}".format(statefile)
        with open('state.json', 'r') as fp:
            ga = pygs.Evolution.from_JSON(fp)
        ga.update_parameters(infile)

    with open('state.json', 'w') as fp:
        fp.write(ga.to_JSON())

    print ga

    print ga.unevaluated_trials

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

    args = parser.parse_args()

    find_groundstate(args.input_file,
                     statefile=args.state)
