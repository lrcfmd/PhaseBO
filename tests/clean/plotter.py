import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm

def get_plot(elements, data, newpoints, dom):
        """
        Plot a contour phase diagram plot, where phase field is colored
        according to degree of instability by interpolation. Currently only
        works for 4-component phase diagrams.

        Returns:
            A matplotlib plot object.
        """

        gridsize = 0.005
        xnew = np.arange(0, 1., gridsize)
        ynew = np.arange(0, 1, gridsize)
 
        f = interpolate.LinearNDInterpolator(data[:, 0:2], data[:, 2])
        znew = np.zeros((len(ynew), len(xnew)))
        for (i, xval) in enumerate(xnew):
            for (j, yval) in enumerate(ynew):
                znew[j, i] = f(xval, yval)

        fig,ax =plt.subplots()
        plt.contourf(xnew, ynew, znew, 15, cmap=cm.gist_heat)
        plt.colorbar()

        plt.scatter(data[:,0], data[:,1], marker='3',c='cyan',lw=1, label=r'Computed points')
        plt.scatter(dom[:,0], dom[:,1], s=20, c='grey', label=r'Possible calculations')
        plt.scatter(newpoints[:,0], newpoints[:,1], s=40, c='lime', label='Suggested next steps') 


        # AXIS
        plt.axis([0,1.,0,1.])
        plt.xlabel(f'{elements[0]}/({elements[0]}+{elements[1]})', fontsize=14)
        plt.ylabel(f'{elements[2]}/({elements[2]}+{elements[3]})', fontsize=14)
        ax.legend(bbox_to_anchor=(0.25, .95),fontsize=10)
        plt.show()
