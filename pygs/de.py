#!/usr/bin/env python

"""
Insert description here.
"""

__author__ = "Alexander Urban"
__date__   = "2014-03-03"

import json
import numpy as np
from pymatgen.io.vaspio import Poscar

class Sublattice(object):

    name  = ''
    sites = []
    occupation = {}

    def __init__(self, name, occupation, site_coords):
        """
        Arguments:
          name (str)           name of the sublattice ('A', 'B', ...)
          occupation (dict)    atom counts for each species, e.g.
                               {'Li' : 4, 'Co' : 4, 'O' : 8}
          site_coords (array)  fractional coordinates of all sites in
                               a 2-d array
        """

        self.name = name
        self.occupation = occupation
        self.sites = np.array(site_coords)

    def __str__(self):
        s  = '\n Instance of Sublattice\n\n'
        s += ' Number of sites : {}\n'.format(self.nsites)
        s += ' Number of atoms : {}\n'.format(self.natoms)
        s += ' Occupation      : ' + str(self.occupation) + '\n'
        return s

    @property
    def nsites(self):
        return len(self.sites)

    @property
    def natoms(self):
        nat = 0
        for element in self.occupation:
            nat += self.occupation[element]
        return nat

class Trial(object):

    decorations = []
    types       = []

    def __init__(self, decorations):
        """
        Arguments:
          decorations (list)  2-d list of atom types for each site
                              (use None for vacancies) and each sublattice

          decorations[i][j] = atom type on j-th site of i-th sublattice
        """

        nsub = len(decorations)
        for i in range(nsub):
            self.decorations.append(np.zeros(len(decorations[i]), dtype=int))
            self.types.append(list(set(decorations[i])))
            for j in range(len(decorations)):
                if decorations[i][j] is None:
                    self.decoration[i][j] = -1
                else:
                    self.decorations[i][j] = self.types.index(decorations[i][j])

    @classmethod
    def from_sublattice(cls, sublattices):
        """
        Arguments:
          sublattices    list of instances of Sublattice
        Returns:
          A Trial with randomized decoration based on the atom counts in
          decoration_dict.
        """

        decorations = []

        for sub in sublattices:
            deco = []
            for element in sub.occupation:
                deco += sub.occupation[element]*[element]
            deco += (sub.nsites - len(deco))*[None]
            decorations.append(deco)

        trial = cls(decorations)
        trial.randomize()

        return trial

    def __str__(self):
        s  = "\n Instance of Trial\n\n"
        s += " Sublattices      : {}\n".format(self.nsublattices)
        s += " Number of sites  : "
        s += (self.nsublattices*"%3d ") % tuple(self.nsites) + "\n"
        s += " Number of atoms  : "
        s += (self.nsublattices*"%3d ") % tuple(self.natoms) + "\n"
        return s

    def __eq__(self, other):
        try:
            eq = True
            for i in range(self.sublattices):
                eq = eq and np.all(self.decorations[i] == other.decorations[i])
        except:
            eq = False
        return eq

    def __ne__(self, other):
        return not self.__eq__(other)

    def __gt__(self, other):
        return NotImplemented
    def __lt__(self, other):
        return NotImplemented
    def __ge__(self, other):
        return NotImplemented
    def __le__(self, other):
        return NotImplemented

    @property
    def nsublattices(self):
        return len(self.decorations)

    @property
    def nsites(self):
        return [len(self.decorations[i]) for i in range(self.nsublattices)]

    @property
    def nvacancies(self):
        return [np.sum(np.where(self.decorations[i]==-1, 1, 0))
                for i in range(self.nsublattices)]

    @property
    def natoms(self):
        return [(self.nsites[i] - self.nvacancies[i])
                for i in range(self.nsublattices)]

    def randomize(self):
        for i in range(self.nsublattices):
            np.random.shuffle(self.decorations[i])

class Evolution(object):

    paramfile = ''

    sitesfile = ''
    avec = np.identity(3)
    sites = []
    site_types = []

    sublattices = []
    size = 0

    trials = []

    def __init__(self, paramfile):
        """
        Arguments:
          paramfile (str)   name of the parameter file in JSON format

        Sample parameter file:

          {
              "size"  : 10,
              "sites" : "sites.vasp",
              "sublattices" : {
                  "A" : { "Li" : 6, "Ni" : 4, "Sb" : 4 },
                  "B" : { "O"  : 20}
              }
          }
        """

        self.paramfile = paramfile

        with open(paramfile, 'r') as fp:
            params = json.load(fp)

        self.size = params['size']
        self.sitesfile = params['sites']

        # pymatgen specific:
        struc = Poscar.from_file(self.sitesfile).structure
        self.avec  = struc.lattice.matrix
        self.sites = np.array(struc.frac_coords)
        self.site_types = np.array([species.symbol for species in struc.species])

        self.sublattices = []
        for name in params['sublattices']:
            occupation  = params['sublattices'][name]
            site_coords = self.sites[self.site_types==name]
            self.sublattices.append(Sublattice(name, occupation, site_coords))

    def __str__(self):
        s  = "\n Instance of Evolution\n\n"
        s += " Population size  : {}\n".format(self.size)
        s += " Sites read from  : {}\n".format(self.sitesfile)
        s += " Num. sublattices : {}\n".format(len(self.sublattices))
        for sub in self.sublattices:
            s += sub.__str__()
        return s

    @property
    def nsublattices(self):
        return len(self.sublattices)

    def create_initial_population(self):
        """
        Randomly generate initial trial vectors
        """
        pass
