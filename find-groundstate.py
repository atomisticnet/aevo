#!/usr/bin/env python

"""
Insert description here.
"""

from __future__ import print_function

__author__ = "Alexander Urban"
__date__   = "2014-03-03"

import argparse

import pygs

#----------------------------------------------------------------------#

def find_groundstate(infile):

    ga = pygs.Evolution(infile)
    print(ga)

#----------------------------------------------------------------------#

if (__name__ == "__main__"):

    parser = argparse.ArgumentParser(
        description     = __doc__+"\n{} {}".format(__date__,__author__),
        formatter_class = argparse.RawDescriptionHelpFormatter )

    parser.add_argument(
        "input_file",
        help    = "Input file in JSON format.")

    args = parser.parse_args()

    find_groundstate(args.input_file)
