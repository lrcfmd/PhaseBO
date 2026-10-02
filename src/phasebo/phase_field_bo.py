import os
import logging
import warnings
import numpy as np
import pandas as pd
import torch
from pymatgen.analysis.phase_diagram import PhaseDiagram
from typing import Optional, Tuple, List

with warnings.catch_warnings():
    # raised by linear_operator (a gpytorch dependency) on torch >= 2.14
    warnings.filterwarnings("ignore", message="`torch.jit.script` is deprecated", category=FutureWarning)
    from botorch.acquisition.logei import qLogExpectedImprovement
    from botorch.fit import fit_gpytorch_mll
    from botorch.generation import MaxPosteriorSampling
    from botorch.models import SingleTaskGP
    from botorch.optim import optimize_acqf_discrete
    from gpytorch.kernels import MaternKernel, ScaleKernel
    from gpytorch.mlls import ExactMarginalLogLikelihood
    from gpytorch.priors import GammaPrior

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
                 acquisition: str = 'qlogei',
                 kernel: str = 'matern',
                 exceptions: Optional[List[str]] = None,
                 allow_negative: bool = False,
                 logger: logging.Logger = None,
                 output_dir: str = '.',
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
        self.acquisition = acquisition
        self.kernel = kernel
        self.next_formulas = next_formulas
        self.exceptions = exceptions
        self.logger = logger or logging.getLogger(__name__)
        self.output_dir = output_dir

        if self.acquisition not in ('qlogei', 'ts'):
            raise ValueError(f'Unsupported acquisition: "{self.acquisition}". Supported: "qlogei", "ts".')
        if self.kernel not in ('matern', 'rbf'):
            raise ValueError(f'Unsupported kernel: "{self.kernel}". Supported: "matern", "rbf".')
        self.logger.info(f"Model: {self.kernel} kernel, {self.acquisition} acquisition, batch {self.batch}")

        self.setBO()
        if self.mode == 'path':
            self.run_path()
        elif self.mode == 'suggest':
            # references are computed but not in X
            self.next = self.suggest_batch(self.X, self.Y, self.unevaluated(self.domain, self.pd_coords))

    def setBO(self) -> None:
        if self.mode == 'path':
            self.logger.info(f"Mode: 'path' with NSEEDS: {self.n_seeds}")
            if self.seeds_type == 'random':
                self.nseeds, self.nseeds_energy = self.get_random_seeds(self.n_seeds, self.exclude)
            elif self.seeds_type == 'segmented':
                self.nseeds, self.nseeds_energy = self.get_seeds_from_segments(self.disect, self.exclude)
            else:
                raise ValueError(f'Unsupported seeds_type: "{self.seeds_type}". Supported: "random", "segmented"')

            X_init = self.nseeds
            Y_init = self.nseeds_energy[:, None]
            self.domain = self.candidates_fc

        elif self.mode == 'suggest':
            # not the references: training on them made path runs find fewer stable compositions
            X_init = self.candidates_fc
            Y_init = self.candidates_energies[:, None]

            if not self.next_formulas:
                self.logger.info("Generating candidate compositions ...")
                self.next_formulas = generate(self.ions, self.formulas, self.exceptions, self.Ntot, self.limits)

            self.domain, self.next_list = self.get_dom_phase()

        elif self.mode == 'generate':
            self.logger.info("Generating candidate compositions, writing to candidates_list.csv")
            self.next_formulas = generate(self.ions, self.formulas, self.exceptions, self.Ntot, self.limits)
            with open(os.path.join(self.output_dir, "candidates_list.csv"), 'a') as cl:
                for f in self.next_formulas:
                    print(f, file=cl)
        else:
            raise ValueError(f'Unsupported mode: "{self.mode}". Supported: "path", "suggest", "generate".')

        if self.mode != 'generate':
            self.X = np.asarray(X_init, dtype=float)
            self.Y = np.asarray(Y_init, dtype=float)

    @staticmethod
    def unevaluated(choices: np.ndarray, X: np.ndarray) -> np.ndarray:
        """Unique rows of choices that are not rows of X."""
        seen = {tuple(x) for x in X}
        choices = np.unique(choices, axis=0)
        return choices[[tuple(c) not in seen for c in choices]]

    def suggest_batch(self, X: np.ndarray, Y: np.ndarray, choices: np.ndarray) -> np.ndarray:
        """Fit a GP to (X, Y) and pick up to self.batch unevaluated points from choices."""
        choices = torch.as_tensor(self.unevaluated(choices, X))
        if len(choices) == 0:
            return np.empty((0, X.shape[1]))

        # BoTorch maximises
        self.model = SingleTaskGP(torch.as_tensor(X), -torch.as_tensor(Y), covar_module=self.covar_module(X.shape[1]))
        fit_gpytorch_mll(ExactMarginalLogLikelihood(self.model.likelihood, self.model))

        q = min(self.batch, len(choices))
        if self.acquisition == 'ts':
            with torch.no_grad():
                X_next = MaxPosteriorSampling(self.model, replacement=False)(choices, num_samples=q)
        else:
            acqf = qLogExpectedImprovement(self.model, best_f=-Y.min())
            X_next, _ = optimize_acqf_discrete(acqf, q=q, choices=choices, unique=True)
        return X_next.detach().numpy()

    def covar_module(self, dims: int) -> Optional[ScaleKernel]:
        """Matérn 5/2 as in the paper (GPyOpt's default), or None for BoTorch's default RBF."""
        if self.kernel == 'rbf':
            return None
        # priors of BoTorch's Matérn default before 0.12
        return ScaleKernel(MaternKernel(nu=2.5, ard_num_dims=dims, lengthscale_prior=GammaPrior(3.0, 6.0)),
                           outputscale_prior=GammaPrior(2.0, 0.15))

    def run_path(self) -> None:
        """Run max_iter batches over the computed phase field."""
        for _ in range(self.iter):
            X_next = self.suggest_batch(self.X, self.Y, self.domain)
            if len(X_next) == 0:
                break
            self.X = np.vstack((self.X, X_next))
            self.Y = np.vstack((self.Y, [[self.f(x)] for x in X_next]))

    def posterior(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Posterior mean and standard deviation (including noise) of the energy above hull at X."""
        with torch.no_grad():
            post = self.model.posterior(torch.as_tensor(X, dtype=torch.float64), observation_noise=True)
            return -post.mean.numpy().ravel(), post.variance.sqrt().numpy().ravel()

    def plot_convergence(self):
        """Plot the energy of each evaluation in the BO path and the best found so far."""
        import matplotlib.pyplot as plt

        n = np.arange(1, len(self.Y) + 1)
        fig, ax = plt.subplots()
        ax.plot(n, self.Y.ravel(), 'o', alpha=0.5, label='Evaluated')
        ax.plot(n, np.minimum.accumulate(self.Y.ravel()), '-', label='Best so far')
        ax.axvline(len(self.nseeds) + 0.5, ls='--', c='grey', label='End of seeds')
        ax.set_xlabel('Evaluation')
        ax.set_ylabel('Energy above hull (meV/atom)')
        ax.legend()
        return plt

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
            self.logger.info(f"{c}  {round(e, 2)}")

        if self.mode == 'path':
            # dicfc names a coordinate's lowest-energy composition, which may not be the seed
            n_seeds = len(self.seeds)
            names = list(self.seeds) + [self.dicfc[self.fcsym(x)][1] for x in self.X[n_seeds:]]
            energies = self.Y.ravel()
            pd.DataFrame({
                'composition': names,
                'energy (meV/atom)': energies.round(2),
                'seed': np.arange(len(names)) < n_seeds,
            }).to_csv(os.path.join(self.output_dir, 'bo_path.csv'), index=False)

            pf = '-'.join(self.elements)
            with open(os.path.join(self.output_dir, f'BO_Path_in_{pf}.txt'), 'a') as f:
                print('Seeds:', file=f)
                print('------', file=f)
                print('Composition     meV/atom above CH', file=f)
                for n, e in zip(names[:n_seeds], energies[:n_seeds]):
                    print(n, round(e, 2), file=f)
                print('\nBO Path:', file=f)
                print('--------', file=f)
                print('Composition     meV/atom above CH', file=f)
                for n, e in zip(names, energies):
                    print(n, round(e, 2), file=f)

        elif self.mode == 'suggest':
            for n in self.next:
                self.logger.info(f"Next: {self.next_list[self.fcsym(n)]}")

    def get_uncertainty(self, mesh=False) -> None:
        """Log standard deviations of surrogate predictions."""
        if mesh:
            bounds = list(zip(self.domain.min(axis=0), self.domain.max(axis=0)))
            X1 = np.linspace(bounds[0][0], bounds[0][1], mesh)
            X2 = np.linspace(bounds[1][0], bounds[1][1], mesh)
            X3 = np.linspace(bounds[2][0], bounds[2][1], mesh)
            x1, x2, x3 = np.meshgrid(X1, X2, X3)
            X = np.hstack((x1.reshape(-1, 1), x2.reshape(-1, 1), x3.reshape(-1, 1)))
            _, std = self.posterior(X)

            self.logger.info(f"Minimum uncertainty in prediction is {round(std.min(), 1)} meV/atom at {X[np.argmin(std)]}")
            self.logger.info(f"Maximum uncertainty in prediction is {round(std.max(), 1)} meV/atom at {X[np.argmax(std)]}")

        else:
            mean, std = self.posterior(self.next_coords)

            un_df = pd.DataFrame({
                'Candidates': self.next_formulas,
                'Posterior mean (meV/atom)': mean.round(1),
                'Posterior std (meV/atom)': std.round(1)
            })
            un_df = un_df.sort_values(['Posterior mean (meV/atom)'])

            path = os.path.join(self.output_dir, 'posterior.csv')
            un_df.to_csv(path, index=False)
            self.logger.info(f"Posterior CSV saved: {path}")
