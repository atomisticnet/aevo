#!/usr/bin/env python

"""
Insert description here.
"""

__author__ = "Alexander Urban"
__date__   = "2014-03-03"

import json
import hashlib
import numpy as np
from pymatgen.io.vaspio import Poscar

class Serializable(object):

    def to_JSON(self):
        def serialize(o):
            if hasattr(o, '__dict__'):
                return o.__dict__
            else:
                try:
                    return o.tolist()
                except:
                    return o
        return json.dumps(self, default=serialize, sort_keys=True, indent=4)

    @classmethod
    def from_JSON(cls, json_string):
        o = cls()
        o.__dict__.update(json.loads(json_string))
        return o

class Sublattice(Serializable):

    def __init__(self, name, occupation, site_index):
        """
        Arguments:
          name (str)           name of the sublattice ('A', 'B', ...)
          occupation (dict)    atom counts for each species, e.g.
                               {'Li' : 4, 'Co' : 4, 'O' : 8}
          site_index (list)    Boolean site index that selects only those sites
                               that belong to this sublattice (ndarray or list)
        """

        self.name = name
        self.occupation = occupation
        self.site_index = np.array(site_index)

    def __str__(self):
        s  = '\n Instance of Sublattice\n\n'
        s += ' Number of sites : {}\n'.format(self.nsites)
        s += ' Number of atoms : {}\n'.format(self.natoms)
        s += ' Occupation      : ' + str(self.occupation) + '\n'
        return s

    @property
    def nsites(self):
        return np.sum(self.site_index)

    @property
    def natoms(self):
        nat = 0
        for element in self.occupation:
            nat += self.occupation[element]
        return nat

class Trial(Serializable):

    def __init__(self, decorations, types=None, fitness=None):
        """
        Arguments:
          decorations (list)  2-d list of atom types for each site
                              (use None for vacancies) and each sublattice
                              i.e. decorations[i][j] = atom type on j-th site
                                                       of i-th sublattice
          types (list)        (optional) Alternatively the atom types can be
                              passed as separate list.  In this case,
                              DECORATIONS is assumed to be a list of integers
                              that refer to the corresponding entry in the
                              TYPES list
          fitness (float)     (optional) Fitness of the trial, if known.
        """

        if types is not None:
            self.types = types
            self.decorations = []
            for deco in decorations:
                self.decorations.append(np.array(deco))
        else:
            self.decorations = []
            self.types       = []
            nsub = len(decorations)
            for i in range(nsub):
                self.decorations.append(np.zeros(len(decorations[i]), dtype=int))
                self.types.append(list(set(decorations[i])))
                for j in range(len(decorations[i])):
                    self.decorations[i][j] = self.types[i].index(decorations[i][j])

        self.fitness = fitness

    @classmethod
    def from_sublattices(cls, sublattices):
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

    def __repr__(self):
        return '<{}(id={})>'.format(self.__class__.__name__, self.id)

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
    def id(self):
        return hashlib.md5(str(self.decorations)).hexdigest()

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

class Evolution(Serializable):

    def __init__(self, avec, sites, site_types):
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

        self.avec = np.array(avec)
        self.sites = np.array(sites)
        self.site_types = np.array(site_types)
        self.sitesfile = None
        self.paramfile = None

        self.sublattices = []
        self.trials = []

    @classmethod
    def from_parameter_file(cls, paramfile):

        with open(paramfile, 'r') as fp:
            params = json.load(fp)

        sitesfile = params['sites']

        # pymatgen specific:
        struc = Poscar.from_file(sitesfile).structure
        avec  = struc.lattice.matrix
        sites = np.array(struc.frac_coords)
        site_types = np.array([species.symbol for species in struc.species])

        evo = cls(avec, sites, site_types)
        evo.sitesfile = sitesfile
        evo.paramfile = paramfile

        for name in params['sublattices']:
            occupation  = params['sublattices'][name]
            site_index = (evo.site_types==name)
            evo.add_sublattice(name, occupation, site_index)

        size = params['size']
        for i in range(size):
            evo.trials.append(Trial.from_sublattices(evo.sublattices))

        return evo

    @classmethod
    def from_JSON(cls, string_or_fp):

        if isinstance(string_or_fp,file):
            entries = json.load(string_or_fp)
        else:
            entries = json.loads(string_or_fp)

        avec = np.array(entries['avec'])
        sites = np.array(entries['sites'])
        site_types = np.array(entries['site_types'])

        evo = cls(avec, sites, site_types)
        evo.sitesfile = entries['sitesfile']
        evo.paramfile = entries['paramfile']

        for sub in entries['sublattices']:
            evo.add_sublattice(**sub)

        for trial in entries['trials']:
            evo.trials.append(Trial(**trial))

        return evo

    def __str__(self):
        s  = "\n Instance of Evolution\n\n"
        if self.paramfile is not None:
            s += " Parameters from  : {}\n".format(self.paramfile)
        if self.sitesfile is not None:
            s += " Sites read from  : {}\n".format(self.sitesfile)
        s += " Population size  : {}\n".format(self.size)
        s += " Num. sublattices : {}\n".format(len(self.sublattices))
        for sub in self.sublattices:
            s += sub.__str__()
        return s

    @property
    def nsublattices(self):
        return len(self.sublattices)

    @property
    def size(self):
        return len(self.trials)

    @property
    def unevaluated_trials(self):
        return [trial for trial in self.trials if trial.fitness is None]

    def update_parameters(self, paramfile):
        """
        Update the algorithm parameters by re-parsing the parameter file.

        Arguments:
          paramfile (str)   path to the parameter file
        """

        self.paramfile = paramfile

        with open(paramfile, 'r') as fp:
            params = json.load(fp)

        sitesfile = params['sites']

        if (self.sitesfile != sitesfile):
            print "Warning: the name of the sites file has changed."
            print "         The new file will be ignored!"

        if (len(params['sublattices']) != self.nsublattices):
            print "Warning: the number of sub-lattices in the input file has changed."
            print "         This update will be ignored!"

        # size = params['size']


    def add_sublattice(self, name, occupation, site_index):
        self.sublattices.append(Sublattice(name, occupation, site_index))

    def write_unevaluated(self, dir='.', format='vasp'):
        """
        Write out structures of trials for whom the fitness has not
        yet been evaluated.

        Arguments:
          dir (str)     path to ooutput directory
          format (str)  atomic structure file format
        """

        pass
