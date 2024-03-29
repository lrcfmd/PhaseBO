import sys
import time
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from numpy.random import seed
from pymatgen.analysis.phase_diagram import PhaseDiagram
from GPyOpt.methods import BayesianOptimization
import matplotlib.pyplot as plt
from matplotlib import cm
from phasebo.phase_field import PhaseField
from phasebo.list_compositions import generate

class PhaseFieldBO(PhaseField):

    def __init__(self,
                 compositions,
                 references, 
                 ions,
                 mode,
                 seeds_type,
                 n_seeds=9,
                 next_formulas=None,
                 exclude_zeros=False,
                 disect=3,
                 Ntot=24,
                 limits=None,
                 max_iter=10, 
                 batch=4,
                 exceptions=None,
                 allow_negative=False):

        super().__init__(compositions, references, ions, exceptions, allow_negative)  
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
        self.setBO()
        if self.mode == 'path':
            self.bo.run_optimization(self.iter, verbosity=False)
        elif self.mode == 'suggest':
            self.next = self.bo.suggest_next_locations()

    def setBO(self):
        if self.mode == 'path':
            print(f"Mode: 'path' with NSEEEDS:{self.n_seeds}")
            if self.seeds_type == 'random':
                self.nseeds, self.nseeds_energy = self.get_random_seeds(self.n_seeds, self.exclude)
            elif self.seeds_type == 'segmented':
                self.nseeds, self.nseeds_energy = self.get_seeds_from_segments(self.disect, self.exclude)
            else:
                raise ValueError(f'Unsupported seeds_type: "{self.seeds_type}". Supported seeds_type are "random" or "segmented"')
            f = self.f
            X_init = self.nseeds
            Y_init = self.nseeds_energy[:,None]
            self.domain = [{'name': 'var_1', 'type': 'bandit', 'domain': self.candidates_fc}]

        elif self.mode == 'suggest':
            f = None
            X_init = self.candidates_fc
            Y_init = self.candidates_energies[:,None]
            if not self.next_formulas:
                print("Generating candidate compositions ...")
                self.next_formulas = generate(self.ions, self.formulas, self.exceptions, self.Ntot, self.limits)
                #for f in self.next_formulas: print (f)

            dom, self.next_list = self.get_dom_phase()
            self.domain = [{'name': 'var_1', 'type': 'bandit', 'domain': dom}]

        elif self.mode == 'generate':
                print("Generating candidate compositions, writing to candidates_list.csv")
                print(self.exceptions)
                self.next_formulas = generate(self.ions, self.formulas, self.exceptions, self.Ntot, self.limits)
                with open("candidates_list.csv",'a') as cl:
                    for f in self.next_formulas:
                        print(f, file=cl)

        else:
           raise ValueError(f'Unsupported mode: "{self.mode}". Supported modes are "path", "suggest" or "generate".') 

        if self.mode != 'generate':
            self.bo = BayesianOptimization(f=f,
                                      domain=self.domain,
                                      X = X_init,
                                      Y = Y_init,
                                      evaluator_type = 'thompson_sampling',
                                      batch_size = self.batch,
                                      de_duplication = True)

    def get_dom_phase(self):
       """ Add generated formulas to the temporary phase field, so their pd_coordinates can be calculated """
       next_entries, next_formulas = self.computed_compositions(self.next_formulas, 100*np.ones(len(self.next_formulas)))
       tmp_pd = PhaseDiagram(self.computed_entries + next_entries)
       self.next_coords = self.get_phase_coordinates(tmp_pd, next_formulas) 
       next_dic = {self.fcsym(f) : name for f, name in zip(self.next_coords, self.next_formulas)}

       return self.next_coords, next_dic 

    def f2f(ex):
        """ TODO: transform atomic fractions (D-1)  to 2D fractions
        add first element (Li) via charge balance  
        """
        return


    def plot_path(self, online=False):
       """ segments considered points to seeds, observed and best.
       if online - plots with annotation on hover """

       observed = [x for x in self.bo.X if all(np.all(self.nseeds != x, axis=1))]
       #observed = [self.f2f(x) for x in observed] 
       color_observed = np.array([self.dicfc[self.fcsym(x)][0] for x in observed])

       best = [self.bo.x_opt]
       for x in self.bo.X:
           if self.dicfc[self.fcsym(x)][0] <= 0:
              best.append(x)
       best = [b for b in best if self.dicfc[self.fcsym(b)][0] <= 0]
       color_best = np.array([self.dicfc[self.fcsym(x)][0] for x in best])

       fig, ax = plt.subplots()

       # transform to 2D:
       #candidates = [self.f2f(x) for x in self.candidates_fc]
       #seeds = [self.f2f(x) for x in self.nseeds]
       
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
        """ works only for 2D square plot with x = C1/(C1+C2) and y = A1/(A1+A2) """
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
       print(f"Writing results to the log-file...")
       f = log
       arg = np.argsort(self.candidates_energies)
       print('All compositions:', file=f)
       print('-----------------', file=f)
       print('Composition     meV/atom above CH', file=f)
       for c, e in zip(np.array(self.candidates)[arg], np.array(self.candidates_energies)[arg]):
            print(c, round(e,2), file=f)
       
       if self.mode == 'path':
           observed = [x for x in self.bo.X] 
           en_observed = np.array([self.dicfc[self.fcsym(x)][0] for x in observed])
           names = [self.dicfc[self.fcsym(x)][1] for x in observed]

           pf = '-'.join(self.elements)
           with open(f'BO_Path_in_{pf}.txt', 'a') as f:
              print('Seeds:', file=f)
              print('------', file=f)
              print('Composition     meV/atom above CH', file=f)
              for s,e in zip(self.nseeds, self.nseeds_energy): 
                  print(self.dicfc[self.fcsym(s)][1], round(e,2), file=f)
              print()
              print('BO Path:', file=f)
              print('--------', file=f) 
              print('Composition     meV/atom above CH', file=f)
              for n, e in zip(names, en_observed):
                  print(n, round(e,2), file=f)

       if self.mode == 'suggest':
           for n in self.next:
               print(n)
               key = self.fcsym(n) 
               #key = f"{n[0]} {n[1]}" # for 2D square
               print(self.next_list[key])
               print("next:", self.next_list[key], file=f)

    def get_uncertainty(self, log, mesh=None):
        """ prints varience of surrogate function """
        model = self.bo.model
        
        if mesh: 
            # for a dense mesh - works for 3 components only 
            bounds = self.bo.acquisition.space.get_bounds()
            X1 = np.linspace(bounds[0][0], bounds[0][1], mesh)
            X2 = np.linspace(bounds[1][0], bounds[1][1], mesh)
            X3 = np.linspace(bounds[2][0], bounds[2][1], mesh)
            x1, x2, x3 = np.meshgrid(X1, X2, X3)
            X = np.hstack((x1.reshape(mesh**3,1),x2.reshape(mesh**3,1),x3.reshape(mesh**3,1))) 
            mean, variance = model.predict(X)
            print(f"Minimum uncertainty in prediction of energy of unexplored compositions is \n \
                   {round(min(variance)[0],1)} meV/atom for coordinates {X[np.argmin(variance)]}", file=log)
            print(f"Maximum uncertainty in prediction of energy of unexplored compositions is \n \
                   {round(max(variance)[0],1)} meV/atom for coordinates {X[np.argmax(variance)]}", file=log)


        else:
            # for all next candidate compositions
            X = self.next_coords
         
            scaler = StandardScaler().fit(self.candidates_energies[:,None])
        
            mean, variance = model.predict(X)
            mean = scaler.inverse_transform(mean)
            variance = scaler.inverse_transform(variance)

            #maxvar = np.argmax(variance)
            #minvar = np.argmin(variance)
            #uncertain = self.next_formulas[maxvar]
            #uncertain_e = round(mean[maxvar][0],1)
            #certain = self.next_formulas[minvar]
            #certain_e = round(mean[minvar][0],1)

            un_df = pd.DataFrame({'Candidates': self.next_formulas, 
                                  'Posterior mean (meV/atom)': [round(i,1) for i in mean.flatten()],
                                  'Variance (meV/atom)': [round(i,1) for i in variance.flatten()]})
            un_df = un_df.sort_values(['Posterior mean (meV/atom)'])

            # --------- logtime
            t = time.localtime()
            timestamp = time.strftime('%b-%d-%Y_%H%M', t)
            un_df.to_csv(f'posterior_{timestamp}.csv', index=None)
