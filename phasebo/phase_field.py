import logging
import numpy as np
import random
from numpy import ndarray
from numpy.random import seed
from typing import List, Dict, Tuple, Optional

from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.core.composition import Composition

class PhaseField:
    """
    Represents a phase field from a list of compositions and their energies.
    Constructs convex hull (phase diagram), extracts fractional coordinates,
    segments the field for seed selection, and provides functions for potential energy surfaces
    to be used in Bayesian optimization.
    """

    def __init__(self,
                 compositions: ndarray,
                 references: ndarray,
                 ions: Dict[str, float],
                 exceptions: Optional[List[str]] = None,
                 allow_negative: bool = True,
                 logger: logging.Logger = None):
        self.logger = logger or logging.getLogger(__name__)
        self.references = list(references[:, 0])
        self.elements = list(ions.keys())
        self.exceptions = exceptions if exceptions else []

        # Will be populated
        self.compositions: List[str] = []
        self.enthalpies: List[float] = []

        # Data containers
        self.pd_coords: ndarray = np.array([])
        self.formulas: List[Composition] = []
        self.computed_entries: List[ComputedEntry] = []
        self.energies: ndarray = np.array([])
        self.dic: Dict[str, List] = {}
        self.dicfc: Dict[str, List] = {}
        self.sections: List[List[str]] = []
        self.seeds: List[str] = []
        self.seeds_energy: List[float] = []
        self.candidates: ndarray = np.array([])
        self.candidates_fc: ndarray = np.array([])
        self.candidates_energies: ndarray = np.array([])

        self.exclude_exceptions(compositions)
        self.compute_convex(allow_negative)
        self.pd_coords = self.get_phase_coordinates(self.pd, self.formulas)
        self.create_dict()
        self.create_dicfc()
        self.get_candidates()

    def exclude_exceptions(self, compositions: ndarray):
        """
        Filter out specified compositions (exceptions) from consideration.
        """
        self.logger.info("Excluding exceptions...")
        for i in range(len(compositions)):
            name = compositions[i, 0].strip()
            if name not in self.exceptions:
                self.compositions.append(name)
                self.enthalpies.append(compositions[i, 1])
            else:
                self.logger.info(f"Excluding from consideration: {name}")

    @staticmethod
    def computed_compositions(compositions: List[str], enthalpies: List[float]) -> Tuple[List[ComputedEntry], List[Composition]]:
        """
        Returns ComputedEntry objects and their corresponding Composition formulas.
        """
        entries, formulas = [], []
        for c, e in zip(compositions, enthalpies):
            ce = ComputedEntry(c, e)
            entries.append(ce)
            formulas.append(ce.composition)
        return entries, formulas

    def compute_convex(self, allow_negative: bool = False):
        """
        Calculates energies above convex hull (meV/atom) for all compositions.
        """
        self.logger.info("Computing energies above convex hull...")
        self.computed_entries, self.formulas = self.computed_compositions(self.compositions, self.enthalpies)
        self.pd = PhaseDiagram(self.computed_entries)
        energies_list = []

        for entry in self.computed_entries:
            hull_energy = 1000 * self.pd.get_e_above_hull(entry)
            if allow_negative and hull_energy == 0:
                try:
                    energies_list.append(1000 * self.pd.get_equilibrium_reaction_energy(entry))
                except Exception as ex:
                    self.logger.info(f"Exception for {entry.composition}: {ex}")
                    energies_list.append(hull_energy)
            else:
                energies_list.append(hull_energy)

        self.energies = np.array(energies_list)

    @staticmethod
    def get_phase_coordinates(pd: PhaseDiagram, formulas: List[Composition]) -> ndarray:
        """
        Get simplex coordinates for compositions in a phase diagram.
        """
        coords = [pd.pd_coords(f) for f in formulas]
        return np.array(coords)

    def get_2D_square_coordinates(self) -> ndarray:
        """
        Compute fractional 2D coordinates for quaternary phase fields.
        """
        var1, var2 = [], []
        for formula in self.formulas:
            c1, c2 = formula[self.elements[0]], formula[self.elements[1]]
            a1, a2 = formula[self.elements[2]], formula[self.elements[3]]
            var1.append(c1 / (c1 + c2) if c1 + c2 else 1)
            var2.append(a1 / (a1 + a2) if a1 + a2 else 1)
        return np.vstack([var1, var2]).T

    def create_dict(self):
        """
        Create dictionary with unique compositions as keys and [energy, fractional coordinates, 2D square coords] as values.
        """
        if len(self.elements) == 4:
            square_coords = self.get_2D_square_coordinates()
        else:
            square_coords = np.zeros((len(self.compositions), 2))

        for c, e, fc, sc in zip(self.compositions, self.energies, self.pd_coords, square_coords):
            if c not in self.dic or e < self.dic[c][0]:
                self.dic[c] = [e, fc, sc]

    @staticmethod
    def fcsym(fc: ndarray) -> str:
        """
        Returns string representation for fractional coordinates.
        """
        return " ".join(map(str, fc.tolist()))

    def create_dicfc(self):
        """
        Create dictionary using fractional coordinate strings as keys and [energy, composition] as values.
        """
        for fc, e, c in zip(self.pd_coords, self.energies, self.compositions):
            self.dicfc[self.fcsym(fc)] = [e, c]

    def get_candidates(self):
        """
        Segregate candidates by excluding references and seeds.
        """
        self.candidates = np.array([c for c in self.compositions if c not in self.references + self.seeds])
        self.candidates_fc = np.array([self.dic[c][1] for c in self.candidates])
        self.candidates_energies = np.array([self.dic[c][0] for c in self.candidates])

    def segment_pf(self, disect: int = 4):
        """
        Segment the 2D phase field into disect x disect sections.
        """
        self.sections = [[[] for _ in range(disect)] for _ in range(disect)]
        rangex = np.linspace(0, 1, disect + 1)
        secx = list(zip(rangex[:-1], rangex[1:]))

        for i, fc in enumerate(self.candidates_fc):
            for x in range(disect):
                if secx[x][0] <= fc[0] < secx[x][1]:
                    for y in range(disect):
                        if secx[y][0] <= fc[1] < secx[y][1]:
                            self.sections[x][y].append(self.candidates[i])

        # Flatten sections for easier access later
        self.sections = [sec for row in self.sections for sec in row]

    def get_seeds_from_segments(self, disect: int = 4, exclude: bool = False) -> Tuple[ndarray, ndarray]:
        """
        Get one random seed from each segment.
        Optionally exclude seeds with zero energies.
        """
        self.segment_pf(disect)

        for sec in self.sections:
            random.shuffle(sec)
            if not sec:
                continue
            if exclude:
                sec = [s for s in sec if self.dic[s][0] != 0]
                if not sec:
                    continue
            self.seeds.append(sec[0])

        self.get_candidates()
        return np.array([self.dic[s][1] for s in self.seeds]), np.array([self.dic[s][0] for s in self.seeds])

    def get_random_seeds(self, n: int, exclude: bool = False) -> Tuple[ndarray, ndarray]:
        """
        Select n random seeds from candidates.
        """
        seeds = self.candidates.tolist()
        random.shuffle(seeds)
        if exclude:
            seeds = [s for s in seeds if self.dic[s][0] != 0]
        self.seeds = seeds[:n]
        self.get_candidates()
        return np.array([self.dic[s][1] for s in self.seeds]), np.array([self.dic[s][0] for s in self.seeds])

    def plot_convex(self):
        """
        Plot interpolated energy above convex hull on the 2D square field.
        """
        import matplotlib.pyplot as plt
        from matplotlib import cm
        from scipy import interpolate

        gridsize = 0.005
        xnew = np.arange(0, 1, gridsize)
        ynew = np.arange(0, 1, gridsize)

        data = self.get_2D_square_coordinates()
        f_interp = interpolate.LinearNDInterpolator(data, self.energies)
        znew = np.array([[f_interp(x, y) for x in xnew] for y in ynew])

        fig, ax = plt.subplots()
        cs = ax.contourf(xnew, ynew, znew, 50, cmap=cm.gist_heat)
        cbar = plt.colorbar(cs, ax=ax)
        cbar.set_label('Energy above hull (meV/atom)', fontsize=14)

        ax.scatter(data[:, 0], data[:, 1], c='lime', marker='3', lw=1.5, label='Computed compositions')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel(f'{self.elements[0]}/({self.elements[0]}+{self.elements[1]})', fontsize=14)
        ax.set_ylabel(f'{self.elements[2]}/({self.elements[2]}+{self.elements[3]})', fontsize=14)
        ax.legend(bbox_to_anchor=(0.5, 1.1), fontsize=10)

        return plt

    def f(self, x: ndarray) -> float:
        """
        Function of energy at fractional coordinate x (discrete).
        """
        matches = np.all(self.pd_coords == x, axis=1)
        return np.dot(np.where(matches, 1, 0), self.energies)
