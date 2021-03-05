import sys
import numpy as np
from numpy.random import seed
from GPyOpt.methods import BayesianOptimization
import matplotlib.pyplot as plt
from matplotlib import cm
from phase_field import PhaseField
import list_compositions

class PhaseFieldBO(PhaseField):

    def __init__(self, compositions, references, ions, mode, seeds_type, exclude_zeros=False, n_seeds=9, disect=3, Ntot=24, max_iter=10, batch=4):
        super().__init__(compositions, references, ions)  
        self.ions = ions
        self.mode = mode
        self.iter = max_iter
        self.domain = []
        self.seeds_type = seeds_type
        self.exclude = exclude_zeros
        self.n_seeds = n_seeds
        self.nseeds = []
        self.nseeds_energy = []
        self.disect = disect
        self.Ntot = Ntot
        self.batch = batch
        self.bo = 0
        self.next = []
        self.next_list = {}
        self.setBO()
        if self.mode == 'path':
            self.bo.run_optimization(self.iter, verbosity=False)
        else:
            self.next = self.bo.suggest_next_locations()

    def setBO(self):
        if self.seeds_type == 'random':
            self.nseeds, self.nseeds_energy = self.get_random_seeds(self.n_seeds, self.exclude)
        elif self.seeds_type == 'segmented':
            self.nseeds, self.nseeds_energy = self.get_seeds_from_segments(self.disect, self.exclude)
        else:
            raise ValueError(f'Unsupported seeds_type: "{self.seeds_type}". Supported seeds_type are "random" or "segmented"')
 
        if self.mode == 'path':
            f = self.f
            X_init = self.nseeds
            Y_init = self.nseeds_energy[:,None]
            self.domain = [{'name': 'var_1', 'type': 'bandit', 'domain': self.candidates_fc}]
        elif self.mode == 'suggest':
            f = None
            X_init = self.candidates_fc
            Y_init = self.candidates_energies[:,None]
            dom, self.next_list = list_compositions.generate(self.ions, self.formulas, self.Ntot)
            self.domain = [{'name': 'var_1', 'type': 'bandit', 'domain': dom}]
        else:
           raise ValueError(f'Unsupported mode: "{self.mode}". Supported modes are "path" or "suggest".') 

        self.bo = BayesianOptimization(f=f,
                                      domain=self.domain,
                                      X = X_init,
                                      Y = Y_init,
                                      evaluator_type = 'thompson_sampling',
                                      batch_size = self.batch,
                                      de_duplication = True)

    def plot_path(self, online=False):
       """ segments considered points to seeds, observed and best.
       if online - plots with annotation on hover """

       observed = [x for x in self.bo.X if all(np.all(self.nseeds != x, axis=1))]
       color_observed = np.array([self.dicfc[self.fcsym(x)][0] for x in observed])
       best = [self.bo.x_opt]
       for x in self.bo.X:
           if self.dicfc[self.fcsym(x)][0] <= 0:
              best.append(x)
       best = [b for b in best if self.dicfc[self.fcsym(b)][0] <= 0]
       color_best = np.array([self.dicfc[self.fcsym(x)][0] for x in best])
       
       fig, ax = plt.subplots()
       sc = ax.scatter(self.candidates_fc[:,0], self.candidates_fc[:,1], s=3, c=self.candidates_energies, cmap=cm.gist_heat, edgecolors='none')
       fig.colorbar(sc)
       ax.scatter(self.nseeds[:,0], self.nseeds[:,1], c=self.nseeds_energy, marker='^', cmap=cm.gist_heat, label='seeds')
       ax.scatter(np.array(observed)[:,0], np.array(observed)[:,1], c=color_observed, cmap=cm.gist_heat, label='checked compositions')
       if len(best) > 0:
           ax.scatter(np.asarray(best)[:,0], np.asarray(best)[:,1], c='red', marker='*', cmap=cm.gist_heat, label='best') 

       if self.seeds_type == 'segmented':
           sections = np.linspace(0, 1, self.disect+1)
           ax.set_xticks(sections, minor=True)     
           ax.set_yticks(sections, minor=True)
           ax.grid(which='minor', linestyle='--', alpha=0.3)

       plt.axis([0,1.,0,1.])
       plt.xlabel(f'{self.elements[0]}/({self.elements[0]}+{self.elements[1]})', fontsize=14)
       plt.ylabel(f'{self.elements[2]}/({self.elements[2]}+{self.elements[3]})', fontsize=14)
       plt.legend()
       plt.show()

       if online:
           import mpld3
           labels = list(self.pf.candidates)
           tooltip = mpld3.plugins.PointLabelTooltip(sc, labels=labels)
           mpld3.plugins.connect(fig, tooltip)
           mpld3.show()

    def plot_suggested(self, online=False):
        fig, ax = plt.subplots()
        sconv = self.plot_convex()
        sconv.scatter(self.next[:,0], self.next[:,1], c='lime', label='suggested next')
        sconv.legend()
        sconv.show()

        if online:
            import mpld3
            labels = list(self.pf.candidates)
            tooltip = mpld3.plugins.PointLabelTooltip(sc, labels=labels)
            mpld3.plugins.connect(fig, tooltip)
            mpld3.show()


    def print_results(self, log):
       arg = np.argsort(self.candidates_energies)
       for c, e in zip(np.array(self.candidates)[arg], np.array(self.candidates_energies)[arg]):
            print(c,e, file=log)
       
       for s,e in zip(self.seeds, self.seeds_energy): 
           print(self.dicfc[self.fcsym(s)][1], e, file=log) 

       if self.mode == 'suggest':
           for n in self.next:
               key = f"{n[0]} {n[1]}"
               print("next:", self.next_list[key], file=log)
