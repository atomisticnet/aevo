#!/usr/bin/env python

"""
Differential Evolution

An implementation of the differential evolution algorithm [1] to
identify near-ground-state atomic orderings.

[1] R. Storn and K. Price, J. Global Optim. 11 (1997) 341–359.

"""

import sys
import os
import json
import hashlib
import numpy as np
import pymatgen as mg
from pymatgen.io.vasp.inputs import Poscar

__author__ = "Alexander Urban"
__date__ = "2014-03-03"


class IncompatibleTrialsException(Exception):
    pass


class PopulationTooSmallException(Exception):
    pass


class TrialsNotEvaluatedException(Exception):
    pass


class NoNewTrialsException(Exception):
    pass


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
        s = '\n Instance of Sublattice\n\n'
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
            self.types = []
            nsub = len(decorations)
            for i in range(nsub):
                self.decorations.append(
                    np.zeros(len(decorations[i]), dtype=int))
                self.types.append(list(set(decorations[i])))
                for j in range(len(decorations[i])):
                    self.decorations[i][j] = self.types[i].index(
                        decorations[i][j])

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
        s = "\n Instance of Trial\n\n"
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
            return self.id == other.id
        except:
            return False

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
    def ntypes(self):
        return [len(self.types[i]) for i in range(self.nsublattices)]

    @property
    def nvacancies(self):
        return [np.sum(np.where(self.decorations[i] == -1, 1, 0))
                for i in range(self.nsublattices)]

    @property
    def natoms(self):
        return [(self.nsites[i] - self.nvacancies[i])
                for i in range(self.nsublattices)]

    def randomize(self):
        for i in range(self.nsublattices):
            np.random.shuffle(self.decorations[i])

    def site_index(self, sublattice):
        """
        Boolean index to select only occupied sites on specified
        sub-lattice.

        Arguments:
          sublattice (int)    ID of the sublattice
        """

        if None in self.types[sublattice]:
            vacID = self.types[sublattice].index(None)
        else:
            vacID = -1
        return np.where(self.decorations[sublattice] != vacID)

    def site_types(self, sublattice):
        """
        Decoration of occupied sites, only.

        Arguments:
          sublattice (int)   ID of the sublattice
        """

        return [self.types[sublattice][i] for i in
                self.decorations[sublattice][self.site_index(sublattice)]]

    def mutate(self, mutation_rate):
        """
        Randomly swap pairs of sites.

        Argument:
          mutation_rate (float)   mutation probability (0 <= x <= 1) for
                                  each site
        """

        for isub in range(self.nsublattices):
            # only consider sublattices with more than one species
            if self.ntypes[isub] <= 1:
                continue
            for s1 in range(self.nsites[isub]):
                r = np.random.random()
                # since we always swap pairs of sites, the actual
                # probability has to be mutation_rate/2
                if r < 0.5*mutation_rate:
                    s2 = s1
                    # pick a second site with different species
                    while (self.decorations[isub][s1] ==
                           self.decorations[isub][s2]):
                        s2 = np.random.randint(self.nsites[isub])

                    # swap the types of the two sites
                    buff = self.decorations[isub][s1]
                    self.decorations[isub][s1] = self.decorations[isub][s2]
                    self.decorations[isub][s2] = buff

    def cross(self, other):
        """
        Cross two trials to form a new one.

        Arguments:
          other (Trial)   the other trial

        Returns:
          a new instance of Trial
        """

        if not (self.nsites == other.nsites and
                self.nsublattices == other.nsublattices and
                self.ntypes == other.ntypes):
            raise IncompatibleTrialsException

        try:
            # create index of type IDs
            other2self = []
            for sl in range(other.nsublattices):
                other2self.append([])
                for t in other.types[sl]:
                    other2self[sl].append(self.types[sl].index(t))
        except:
            raise IncompatibleTrialsException

        """
        Crossing algorithm

        trial 1:  [ 1 | 0 | 2 | 1 | 1 | 0 | 2 ]
                    ^       ^   ^       ^
        trial 2:  [ 0 | 1 | 2 | 1 | 2 | 1 | 0 ]
                        ^           ^       ^
        combined: [ 1 | 0 | 2 | 1 | 2 | 0 | 0 ]
        --> one 0 too many, one 1 too few
            replace one 0 with a 1 to fix stoichiometry

        result:   [ 1 | 0 | 2 | 1 | 2 | 0 | 1 ]
        """

        # (1) start with copy of trial 1
        offspring = Trial(self.decorations, types=self.types)

        # (2) replace roughly half of occupations with trial 2
        for sl in range(self.nsublattices):
            for i in range(self.nsites[sl]):
                r = np.random.random()
                if (r > 0.5):  # = 50% probability
                    site_type = other2self[sl][other.decorations[sl][i]]
                    offspring.decorations[sl][i] = site_type

        # (3) fix the stoichiometry
        for sl in range(self.nsublattices):
            diff = []
            required = []
            for i in range(self.ntypes[sl]):
                n1 = np.sum(self.decorations[sl] == i)
                n2 = np.sum(offspring.decorations[sl] == i)
                diff.append(n2-n1)
                if (n1 - n2 > 0):
                    required += (n1-n2)*[i]
            for i in range(self.ntypes[sl]):
                for j in range(diff[i]):
                    sites_i = np.arange(offspring.nsites[sl]
                                        )[offspring.decorations[sl] == i]
                    r1 = np.random.randint(len(sites_i))
                    sel = sites_i[r1]
                    r2 = np.random.randint(len(required))
                    offspring.decorations[sl][sel] = required.pop(r2)

        return offspring


class Evolution(Serializable):

    def __init__(self, avec, sites, site_types, atom_types, size,
                 generation=0):
        """
        Arguments:
          avec (2d array/list)   lattice vectors
          sites (2d array/list)  fractional site coordinates
          site_types (list)      list of sublattice types
          atom_types (list)      list of atomic species (determines order)
          size (int)             number of trials per generation
          generation (int)       generation of the evolutionary algorithm
        """

        self.avec = np.array(avec)
        self.sites = np.array(sites)
        self.site_types = np.array(site_types)
        self.atom_types = atom_types
        self.size = size
        self.generation = generation
        self.sitesfile = None
        self.paramfile = None

        self.sublattices = []
        self.trials = []
        self.fitness_history = {}

    @classmethod
    def from_parameter_file(cls, paramfile):
        """
        Arguments:
          paramfile (str)   path to parameter file

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

        with open(paramfile, 'r') as fp:
            params = json.load(fp)

        sitesfile = params['sites']
        atom_types = params['atom_types']

        # pymatgen specific
        struc = Poscar.from_file(sitesfile).structure
        avec = struc.lattice.matrix
        sites = np.array(struc.frac_coords)
        site_types = np.array([species.symbol for species in struc.species])

        size = params['size']

        evo = cls(avec, sites, site_types, atom_types, size)
        evo.sitesfile = sitesfile
        evo.paramfile = paramfile

        for name in params['sublattices']:
            occupation = params['sublattices'][name]
            site_index = (evo.site_types == name)
            evo.add_sublattice(name, occupation, site_index)

        if 'initial_population' in params:
            for filename in params['initial_population']:
                print " Adding trial structure from file: {}".format(filename)
                # pymatgen specific
                struc = Poscar.from_file(filename).structure
                avec = struc.lattice.matrix
                sites = np.array(struc.frac_coords)
                types = np.array([species.symbol for species in struc.species])
                evo.add_trial(avec, sites, types)
                if (evo.ntrials >= size):
                    break
            print

        for i in range(size-evo.ntrials):
            evo.trials.append(Trial.from_sublattices(evo.sublattices))

        evo.set_algorithm_parameters(params)

        return evo

    @classmethod
    def from_JSON(cls, string_or_fp):

        if isinstance(string_or_fp, file):
            entries = json.load(string_or_fp)
        else:
            entries = json.loads(string_or_fp)

        avec = np.array(entries['avec'])
        sites = np.array(entries['sites'])
        site_types = np.array(entries['site_types'])
        atom_types = entries['atom_types']
        size = entries['size']
        generation = entries['generation']

        evo = cls(avec, sites, site_types, atom_types,
                  size, generation=generation)
        evo.sitesfile = entries['sitesfile']
        evo.paramfile = entries['paramfile']
        evo.fitness_history = entries['fitness_history']

        for sub in entries['sublattices']:
            evo.add_sublattice(**sub)

        for trial in entries['trials']:
            evo.trials.append(Trial(**trial))

        return evo

    def __str__(self):
        s = "\n Instance of Evolution\n\n"
        if self.paramfile is not None:
            s += " Parameters from  : {}\n".format(self.paramfile)
        if self.sitesfile is not None:
            s += " Sites read from  : {}\n".format(self.sitesfile)
        s += " Population size  : {}\n".format(self.ntrials)
        s += " Generation size  : {}\n".format(self.size)
        s += " Num. sublattices : {}\n".format(len(self.sublattices))
        for sub in self.sublattices:
            s += sub.__str__()
        return s

    @property
    def nsublattices(self):
        return len(self.sublattices)

    @property
    def nsites(self):
        return len(self.sites)

    @property
    def ntrials(self):
        return len(self.trials)

    @property
    def best_trial(self):
        fitness = [trial.fitness for trial in self.evaluated_trials]
        idx = np.argsort(fitness)
        return self.trials[idx[0]]

    @property
    def unevaluated_trials(self):
        return [trial for trial in self.trials if trial.fitness is None]

    @property
    def evaluated_trials(self):
        return [trial for trial in self.trials if trial.fitness is not None]

    @property
    def average_fitness(self):
        fitness = [trial.fitness for trial in self.evaluated_trials]
        return np.sum(fitness)/len(fitness)

    @property
    def min_fitness(self):
        fitness = [trial.fitness for trial in self.evaluated_trials]
        return np.min(fitness)

    @property
    def max_fitness(self):
        fitness = [trial.fitness for trial in self.evaluated_trials]
        return np.max(fitness)

    def sublattice_of_site(self, isite):
        """
        Return sublattice ID of site ISITE.
        """
        return [i for i in range(self.nsublattices)
                if self.sublattice[i].name == self.site_types[isite]][0]

    def top_N(self, N):
        """
        Return list of tuples (ID, fitness) of the N best trials found so far.
        """
        return [(t, self.fitness_history[t]) for t in
                sorted(self.fitness_history,
                       key=self.fitness_history.get)[0:N]]

    def update_parameters(self, paramfile):
        """
        Update the algorithm parameters by re-parsing the parameter file.

        Arguments:
          paramfile (str)  path to the parameter file
        """

        with open(paramfile, 'r') as fp:
            params = json.load(fp)

        sitesfile = params['sites']

        if (self.sitesfile != sitesfile):
            print("Warning: the name of the sites file has changed.")
            print("         The new file will be ignored!")

        if (len(params['sublattices']) != self.nsublattices):
            print("Warning: the number of sub-lattices in the input "
                  "file has changed.")
            print("         This update will be ignored!")

        self.paramfile = paramfile

        self.set_algorithm_parameters(params)

    def set_algorithm_parameters(self, params):
        """
        Set parameters of the evolutionary algorithm.

        Arguments:
          params (dict)  dictionary with parameters, e.g., from
                         JSON file
        """

        self.size = params['size']
        self.keep = params['keep']
        self.mutate = params['mutate']

    def add_sublattice(self, name, occupation, site_index):
        self.sublattices.append(Sublattice(name, occupation, site_index))

    def add_trial(self, avec, sites, types):
        """
        Add a trial structure to the current population.

        Arguments:
          avec (2-d array)    lattice vectors
          sites (2d array)    site coordinates
          types (list)        list of types of all sites (species)
        """

        # For now the lattice vectors are ignored, and it is assumed
        # that a compatible lattice is used.  No check here, to allow
        # for small cell relaxations.

        decorations = []
        for subl in self.sublattices:
            decorations.append(subl.nsites*[None])

        for i1 in range(len(sites)):
            s1 = sites[i1]
            d_min = np.linalg.norm(self.avec[0])
            d_min += np.linalg.norm(self.avec[1])
            d_min += np.linalg.norm(self.avec[2])
            isub_min = -1
            i2_min = -1
            for isub in range(self.nsublattices):
                subl = self.sublattices[isub]
                for i2 in range(subl.nsites):
                    s2 = self.sites[subl.site_index][i2]
                    d = np.linalg.norm(s2 - s1)
                    if (d < d_min):
                        d_min = d
                        isub_min = isub
                        i2_min = i2
            if (isub_min >= 0) and (i2_min >= 0):
                if decorations[isub_min][i2_min] is None:
                    decorations[isub_min][i2_min] = types[i1]
                else:
                    sys.stderr.write("Error: overlapping sites detected.\n")
                    raise IncompatibleTrialsException
            else:
                sys.stderr.write(
                    "Error: incompatible trial structure skipped.\n")
                raise IncompatibleTrialsException

        self.trials.append(Trial(decorations))

    def trial_coords(self, trial, sort=True):
        """
        Get site coordinates and decoration of given trial.

        Arguments:
          trial (Trial)

        Returns:
          tuple (coords, types) with ndarrays.
        """

        coords = []
        types = []
        for i in range(self.nsublattices):
            idx = self.sublattices[i].site_index
            coo = self.sites[idx]
            coo = coo[trial.site_index(i)]
            coords += list(coo)
            types += list(trial.site_types(i))
        coords = np.array(coords)
        types = np.array(types)

        if sort:
            coords_sorted = []
            types_sorted = []
            for t in self.atom_types:
                idx = (types == t)
                coords_sorted.extend(coords[idx])
                types_sorted.extend(types[idx])
            coords = np.array(coords_sorted)
            types = np.array(types_sorted)

        return (coords, types)

    def write_trial(self, trial, filename=None, directory='.', frmt='vasp'):
        """
        Save atomic structure of a trial.

        Arguments:
          trial (Trial)
          filename (str)  name of the output file
          directory (str) path to ooutput directory
          frmt (str)      atomic structure file format
        """

        if frmt != 'vasp':
            sys.stderr.write("Error: only VASP format implemented\n")
            raise NotImplementedError

        (coords, types) = self.trial_coords(trial)

        if filename is None:
            filename = directory + os.sep + 'POSCAR-{}.vasp'.format(trial.id)
        print "   writing file: {}".format(filename)

        # pymatgen specific
        s = mg.Structure(lattice=self.avec, species=types, coords=coords)
        p = Poscar(structure=s)
        p.write_file(filename)

    def write_unevaluated(self, directory=None, frmt='vasp'):
        """
        Write out structures of trials for whom the fitness has not
        yet been evaluated.

        Arguments:
          directory (str)   path to output directory
          frmt (str)        atomic structure file format
        """

        print " Writing non-evaluated trial structures to files."
        print

        if directory is None:
            directory = 'generation{:05d}'.format(self.generation)

        # fixme: (1) could be file, not dir (2) racing condition
        if not os.path.exists(directory):
            os.makedirs(directory)

        for trial in self.unevaluated_trials:
            self.write_trial(trial, directory=directory, frmt=frmt)

        print

    def write_current_best(self, filename, directory='.', frmt='vasp'):
        """
        Write out structures of the current best trial.

        Arguments:
          filename (str)    path to the output file
          directory (str)   path to output directory
          frmt (str)        atomic structure file format
        """

        print " Saving current best trial structures to file."
        print

        # fixme: (1) could be file, not dir (2) racing condition
        if not os.path.exists(directory):
            os.makedirs(directory)

        self.write_trial(self.best_trial, filename=filename,
                         directory=directory, frmt=frmt)

        print

    def read_fitness(self, directory=None):
        """
        Read fitness values of current trials from files.

        Arguments:
          directory (str)   path to directory with FITNESS files
        """

        print " Reading fitness values."
        print

        if directory is None:
            directory = 'generation{:05d}'.format(self.generation)

        for trial in self.unevaluated_trials:
            filename = directory + os.sep + 'FITNESS.' + trial.id
            if os.path.exists(filename):
                print "   reading file: {}".format(filename)
                with open(filename, 'r') as fp:
                    trial.fitness = float(fp.readline())
                self.fitness_history[trial.id] = trial.fitness
            else:
                print "   file not found: {}".format(filename)

        print

    def select(self, keep=None):
        """
        Select trials for the next generation.

        Arguments:
          n (int)       Number of trials to be selected.
          keep (int)    Number of current best trials to keep in
                        the next generation.  Further trials will
                        be selected probabilistically.
        """

        if keep is None:
            keep = self.keep

        # check if population is large enough
        if (self.ntrials < self.size):
            raise PopulationTooSmallException

        # check if all trials have been evaluated
        if (len(self.unevaluated_trials) > 0):
            raise TrialsNotEvaluatedException

        self.generation += 1

        if (self.ntrials == self.size):
            # no selection required
            return

        keep = min(keep, self.size)
        new_population = []

        fitness = [trial.fitness for trial in self.trials]
        rank = np.argsort(fitness)
        for i in range(keep):
            new_population.append(self.trials[rank[i]])
        for i in range(keep):
            self.trials.remove(new_population[i])

        # Roulette Wheel selection of further trials
        while(self.size > len(new_population)):
            # (1) score each remaining trial proportionally to its fitness
            fitness = [trial.fitness for trial in self.trials]
            rank = np.argsort(fitness)
            fit_min = fitness[rank[0]]   # best
            fit_max = fitness[rank[-1]]  # worst
            score = [fit_max - fit_min]
            for i in range(1, self.ntrials):
                score.append(fit_max - self.trials[i].fitness + score[i-1])
            # (2) select trials
            r = np.random.random()*score[-1]
            i = 0
            while(score[i] < r):
                i += 1
            new_population.append(self.trials.pop(i))

        # delete un-selected trials
        for trial in self.trials:
            self.trials.remove(trial)

        self.trials = new_population

    def mate(self, mutate=None):
        """
        Combine the trials of the current generation to generate the
        next generation of trials.

        Arguments:
          mutate (float in [0.0,1.0])   mutation probability
        """

        if mutate is None:
            mutate = self.mutate

        mutate = min(self.nsites, mutate)

        # check if population is large enough
        if (self.ntrials < self.size):
            raise PopulationTooSmallException

        # check if all trials have been evaluated
        if (len(self.unevaluated_trials) > 0):
            sys.stderr.write("Error: not all trials evaluated."
                             "  Can not create new population.")
            raise TrialsNotEvaluatedException

        new_generation = []
        ntry = 0
        ntry_max = 100
        while (len(new_generation) < self.size):
            # (1) select 2 trials from current generation
            t1 = np.random.randint(self.ntrials)
            t2 = np.random.randint(self.ntrials - 1)
            if (t2 >= t1):
                t2 += 1
            # (2) cross these two
            trial = self.trials[t1].cross(self.trials[t2])
            # (3) introduce mutations
            trial.mutate(mutate)
            # (4) add to new generation, if the trial is truely new
            if ((trial.id not in self.fitness_history) and
                    (trial.id not in [t.id for t in new_generation])):
                new_generation.append(trial)
                ntry = 0
            else:
                ntry += 1
            if (ntry > ntry_max):
                raise NoNewTrialsException

        self.trials.extend(new_generation)
