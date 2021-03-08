import numpy as np
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.entries.computed_entries import ComputedEntry
import random
from numpy.random import seed

class PhaseField:
    """ 
    For the list of compositions and their energies
    determines fractional positions and a phase diagram (convex hull)
    segments the field into dxd sections and select seeds in them,
    defines a function of Potential energy surface for Bayesian Optimisation"""

    def __init__(self, compositions, references, ions):
        self.compositions = compositions[:,0]
        self.enthalpies = compositions[:,1]
        self.references = list(references[:,0])
        self.elements = list(ions.keys())
        self.pd_coords = []
        self.dic = {}
        self.dicfc = {}
        self.sections = []
        self.seeds = []
        self.seeds_energy = []
        # for each instance of a phase field:
        self.compute_convex()
        self.pd_coords = self.get_phase_coordinates(self.pd, self.formulas)
        self.create_dict()
        self.create_dicfc()
        self.get_candidates()

    def f(self, x):
         """ function of ehull energies for all fractional coordinates """
         return np.dot(np.where(np.all(self.pd_coords == x, axis=1), 1, 0), self.energies) 

    @staticmethod
    def computed_compositions(compositions, enthalpies):
        """ Returns a list of ComputedEntry pymatgen objects and a list of ordered by element formulas """
        structures = []
        formulas = []
        for c, e in zip(compositions, enthalpies):
            ce = ComputedEntry(f"{c}", e )
            structures.append(ce)
            formulas.append(ce.composition)
        return structures, formulas

    def compute_convex(self):
        """ Calculates energies above the convex hull (meV/atom)
        for all compositions, including references. 
        Rewrites energies"""
        print ("Computing energies above convex hull...") 
        self.computed_entries, self.formulas = self.computed_compositions(self.compositions, self.enthalpies)
        self.pd = PhaseDiagram(self.computed_entries)
        energies_list = []
        for entry in self.computed_entries:
            energies_list.append(1000*self.pd.get_e_above_hull(entry))
        self.energies = np.array(energies_list)

    @staticmethod
    def get_phase_coordinates(pd, formulas):
        """ Computes coordinates in a simplex of a phase field """
        coords = []
        for formula in formulas:
            coords.append(pd.pd_coords(formula))
        return np.array(coords)
   
    def get_2D_square_coordinates(self):
        """ Computes fractional coordinates of compositions on a 2D(quaternary) phase field""" 
        var1, var2 = [], []
        for formula in self.formulas:
            c1 = formula[self.elements[0]]
            c2 = formula[self.elements[1]]
            a1 = formula[self.elements[2]]
            a2 = formula[self.elements[3]]
            if c1 + c2 != 0:
                var1.append(c1 / (c1 + c2))
            else:
                var1.append(1)
            if a1 + a2 != 0:
                var2.append(a1 / (a1 + a2))
            else:
                var2.append(1)
        return np.array([np.array(var1), np.array(var2)]).T 

    def create_dict(self): 
        """ get a dictionary of unique compositions with lowest energy values
            and fractional coordinates"""
        # for ph. f. representation on a square 
        if len(self.elements) == 4:
            square_coordinates = self.get_2D_square_coordinates()
        else:
            square_coordinates = np.zeros(len(self.compositions))

        for c, e, f, s in zip(self.compositions, self.energies, self.pd_coords, square_coordinates):
            if c not in self.dic: self.dic[c] = [e, f, s]
            elif e < self.dic[c][0]: self.dic[c] = [e, f, s]

    @staticmethod
    def fcsym(fc):
        """ return a string for fractional coordinates """
        return ' '.join(map(str, list(fc)))

    def create_dicfc(self):
        """ dictionary with fractional coordinates for the keys """
        for f, e, c in zip(self.pd_coords, self.energies, self.compositions):
            self.dicfc[self.fcsym(f)] = [e, c]

    def get_candidates(self):
        """ segregate candidates from references and seeds
        and get fractional coordinates of the candidate compositions"""
        self.candidates = np.array([c for c in self.compositions if c not in self.references + list(self.seeds)])
        self.candidates_fc = np.array([self.dic[c][1] for c in self.candidates])
        self.candidates_energies = np.array([self.dic[c][0] for c in self.candidates])

    def segment_pf(self, disect=4):             
        """ segment quaternary phase field (2d) by disect x disect sections. """ 
        sections = [[[] for i in range(disect)] for j in range(disect)]
        rangex = np.linspace(0, 1, disect + 1)
        secx = list(zip(rangex, rangex[1:]))

        for i, c in enumerate(self.candidates_fc):
            for x in range(disect):
                for y in range(disect):
                    if secx[x][0] <= c[0] < secx[x][1]:
                        if secx[y][0] <= c[1] < secx[y][1]:
                            sections[x][y].append(self.candidates[i])
     
        for sec in sections:
            self.sections += sec
   
    def get_seeds_from_segments(self, disect=4, exclude=False):
        """ get random seed from each segment. they will be excluded from candidates.
        return fractional coordinates for the seeds and their energies"""
        self.segment_pf(disect)

        for sec in self.sections:
            random.shuffle(sec)
            if exclude:
                while not self.dic[sec[0]][0]: random.shuffle(sec)
            self.seeds.append(sec[0])

        self.get_candidates()
        return np.array([self.dic[s][1] for s in self.seeds]), np.array([self.dic[s][0] for s in self.seeds]), 

    def get_random_seeds(self, n, exclude=False):
        """ choose seeds randomly. can be called instead of get_seeds_from_segments,
        but not after in the same instance as that excludes more seeds from candidates.
        if exclude==True, excludes seeds that have zero energies """
        seeds = self.candidates
        random.shuffle(seeds)
        if exclude:
            while not all([self.dic[s][0] for s in seeds[:n]]): random.shuffle(seeds)
        self.seeds = seeds[:n]
        self.get_candidates()
        return np.array([self.dic[s][1] for s in self.seeds]), np.array([self.dic[s][0] for s in self.seeds])

    def plot_convex(self):
        """ plot interpolated ehull energy """

        from scipy import interpolate
        from matplotlib import cm
        import matplotlib.pyplot as plt

        gridsize = 0.005
        xnew = np.arange(0, 1., gridsize)
        ynew = np.arange(0, 1, gridsize)
        data = self.get_2D_square_coordinates()
        f = interpolate.LinearNDInterpolator(data, self.energies)
        znew = np.zeros((len(ynew), len(xnew)))

        for (i, xval) in enumerate(xnew):
            for (j, yval) in enumerate(ynew):
                znew[j, i] = f(xval, yval) 

        fig,ax =plt.subplots()
        plt.contourf(xnew, ynew, znew, 50, cmap=cm.gist_heat)
        plt.colorbar()
        plt.scatter(data[:,0], data[:,1], c='cyan', marker='3', lw=0.5, label='Computed compositions')
        plt.axis([0,1.,0,1.])
        plt.xlabel(f'{self.elements[0]}/({self.elements[0]}+{self.elements[1]})', fontsize=14)
        plt.ylabel(f'{self.elements[2]}/({self.elements[2]}+{self.elements[3]})', fontsize=14)
        ax.legend(bbox_to_anchor=(0.5, 1.0),fontsize=10)
        return plt

if __name__ == "__main__":
    import pandas

    ifile = 'LiSnSCl_700eV.csv'           # File with compositions and total enthalpies
    ions = {'Li':1,'Sn':4,'S':-2,'Cl':-1} # Ions and oxidation states 
    df = pandas.read_csv(ifile, header=0)     #
    compositions = df.values              # Select candidate and reported compositions
    references = df.values[195:]          #

    pf = PhaseField(compositions, references, ions)
    coords = pf.pd_coords
    for i, c in enumerate(coords):
        print(pf.formulas[i], c) 

    plot = pf.plot_convex() 
    plot.show()
