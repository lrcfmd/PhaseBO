import sys
import time
import logging
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from numpy.random import seed
from pymatgen.analysis.phase_diagram import PhaseDiagram
from GPyOpt.methods import BayesianOptimization
import matplotlib.pyplot as plt
from matplotlib import cm
from typing import Optional, Tuple, List

from phasebo.phase_field import PhaseField
from phasebo.list_compositions import generate

class PhaseFieldBO(PhaseField):

    def __init__(self,
                 compositions,
                 references,
                 ions,
                 mode: str = 'suggest',
                 seeds_type: str = 'random',
                 n_seeds: int = 9,
                 next_formulas: Optional[List[str]] = None,
                 exclude_zeros: bool = False,
                 disect: int = 3,
                 Ntot: int = 24,
                 limits: Optional[List[float]] = None,
                 max_iter: int = 10,
                 batch: int = 4,
                 exceptions: Optional[List[str]] = None,
                 allow_negative: bool = False,
                 logger: logging.Logger = None,
                 ) -> None:

        super().__init__(compositions, references, ions, exceptions, allow_negative, logger)
        self.ions = ions
        self.mode = mode
        self.iter = max_iter
        self.seeds_type = seeds_type
        self.n_seeds = n_seeds
        self.exclude = exclude_zeros
        self.disect = disect
        self.Ntot = Ntot
        self.limits = limits
        self.batch = batch
        self.next_formulas = next_formulas
        self.exceptions = exceptions
        self.logger = logger or logging.getLogger(__name__)

        self.setBO()
        if self.mode == 'path':
            self.bo.run_optimization(self.iter, verbosity=False)
        elif self.mode == 'suggest':
            self.next = self.bo.suggest_next_locations()

    def setBO(self) -> None:
        if self.mode == 'path':
            self.logger.info(f"Mode: 'path' with NSEEDS: {self.n_seeds}")
            if self.seeds_type == 'random':
                self.nseeds, self.nseeds_energy = self.get_random_seeds(self.n_seeds, self.exclude)
            elif self.seeds_type == 'segmented':
                self.nseeds, self.nseeds_energy = self.get_seeds_from_segments(self.disect, self.exclude)
            else:
                raise ValueError(f'Unsupported seeds_type: "{self.seeds_type}". Supported: "random", "segmented"')

            f = self.f
            X_init = self.nseeds
            Y_init = self.nseeds_energy[:, None]
            self.domain = [{'name': 'var_1', 'type': 'bandit', 'domain': self.candidates_fc}]

        elif self.mode == 'suggest':
            f = None
            X_init = self.candidates_fc
            Y_init = self.candidates_energies[:, None]

            if not self.next_formulas:
                self.logger.info("Generating candidate compositions ...")
                self.next_formulas = generate(self.ions, self.formulas, self.exceptions, self.Ntot, self.limits)

            dom, self.next_list = self.get_dom_phase()
            self.domain = [{'name': 'var_1', 'type': 'bandit', 'domain': dom}]

        elif self.mode == 'generate':
            self.logger.info("Generating candidate compositions, writing to candidates_list.csv")
            self.next_formulas = generate(self.ions, self.formulas, self.exceptions, self.Ntot, self.limits)
            with open("candidates_list.csv", 'a') as cl:
                for f in self.next_formulas:
                    print(f, file=cl)
        else:
            raise ValueError(f'Unsupported mode: "{self.mode}". Supported: "path", "suggest", "generate".')

        if self.mode != 'generate':
            self.bo = BayesianOptimization(f=f,
                                           domain=self.domain,
                                           X=X_init,
                                           Y=Y_init,
                                           evaluator_type='thompson_sampling',
                                           batch_size=self.batch,
                                           de_duplication=True)

    def get_dom_phase(self) -> Tuple[np.ndarray, dict]:
        """Add generated formulas to phase field to compute coordinates."""
        next_entries, next_formulas = self.computed_compositions(self.next_formulas, 100 * np.ones(len(self.next_formulas)))
        tmp_pd = PhaseDiagram(self.computed_entries + next_entries)
        self.next_coords = self.get_phase_coordinates(tmp_pd, next_formulas)
        next_dic = {self.fcsym(f): name for f, name in zip(self.next_coords, self.next_formulas)}
        return self.next_coords, next_dic

    def print_results(self) -> None:
        self.logger.info("Writing results to log file...")
        arg = np.argsort(self.candidates_energies)
        self.logger.info('All compositions:')
        self.logger.info('-----------------')
        self.logger.info('Composition     meV/atom above CH')
        for c, e in zip(np.array(self.candidates)[arg], np.array(self.candidates_energies)[arg]):
            self.logger.info(c, round(e, 2))

        if self.mode == 'path':
            observed = self.bo.X
            en_observed = np.array([self.dicfc[self.fcsym(x)][0] for x in observed])
            names = [self.dicfc[self.fcsym(x)][1] for x in observed]

            pf = '-'.join(self.elements)
            with open(f'BO_Path_in_{pf}.txt', 'a') as f:
                print('Seeds:', file=f)
                print('------', file=f)
                print('Composition     meV/atom above CH', file=f)
                for s, e in zip(self.nseeds, self.nseeds_energy):
                    print(self.dicfc[self.fcsym(s)][1], round(e, 2), file=f)
                print('\nBO Path:', file=f)
                print('--------', file=f)
                print('Composition     meV/atom above CH', file=f)
                for n, e in zip(names, en_observed):
                    print(n, round(e, 2), file=f)

        elif self.mode == 'suggest':
            for n in self.next:
                self.logger.info(f"Next: {self.next_list[self.fcsym(n)]}")

    def get_uncertainty(self, mesh=False) -> None:
        """Log variances of surrogate predictions."""
        model = self.bo.model

        if mesh:
            bounds = self.bo.acquisition.space.get_bounds()
            X1 = np.linspace(bounds[0][0], bounds[0][1], mesh)
            X2 = np.linspace(bounds[1][0], bounds[1][1], mesh)
            X3 = np.linspace(bounds[2][0], bounds[2][1], mesh)
            x1, x2, x3 = np.meshgrid(X1, X2, X3)
            X = np.hstack((x1.reshape(-1, 1), x2.reshape(-1, 1), x3.reshape(-1, 1)))
            _, variance = model.predict(X)

            self.logger.info(f"Min variance: {round(min(variance)[0], 1)} meV/atom at {X[np.argmin(variance)]}")
            self.logger.info(f"Max variance: {round(max(variance)[0], 1)} meV/atom at {X[np.argmax(variance)]}")
            self.logger.info(f"Minimum uncertainty in prediction is {round(min(variance)[0], 1)} meV/atom at {X[np.argmin(variance)]}")
            self.logger.info(f"Maximum uncertainty in prediction is {round(max(variance)[0], 1)} meV/atom at {X[np.argmax(variance)]}")

        else:
            X = self.next_coords
            scaler = StandardScaler().fit(self.candidates_energies[:, None])
            mean, variance = model.predict(X)
            mean = scaler.inverse_transform(mean)
            variance = scaler.inverse_transform(variance)

            un_df = pd.DataFrame({
                'Candidates': self.next_formulas,
                'Posterior mean (meV/atom)': [round(i, 1) for i in mean.flatten()],
                'Variance (meV/atom)': [round(i, 1) for i in variance.flatten()]
            })
            un_df = un_df.sort_values(['Posterior mean (meV/atom)'])

            timestamp = time.strftime('%b-%d-%Y_%H%M', time.localtime())
            un_df.to_csv(f'posterior_{timestamp}.csv', index=False)
            self.logger.info(f"Posterior CSV saved: posterior_{timestamp}.csv")
