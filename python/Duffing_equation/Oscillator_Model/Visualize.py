#
# Visualization class
#

import numpy as np
import matplotlib.pyplot as plt

class Visualize:
    """
    A class for Visualizing the oscillator equation for ion channels
    Parameters
    ----------
    model: Model class
    solution: Solution class
        Solution(model, tmax=1, nsteps=1000)
    """
    def __init__(self, model, solution):
        self.model = model
        self.solution = solution
        
    def plot(self, ratio = 0.9):
        """
        Generate a figure of the drug quantity per
        compartment over time for the corresponding model
        :param separate: set to True if you want 1 plot per compartment
        :returns: matplotlib figure
        """
        solution = self.solution.sol
        subset = int(self.solution.nsteps * ratio)
        fig = plt.figure()
        plt.ylabel('current',fontsize=18)
        plt.xlabel('time',fontsize=18)
        plt.plot(solution.t[subset:], solution.y[0, subset:])
        fig.tight_layout()
        return fig
    
    def plot_poincare(self, ratio=0, name=0, pixel=1):
        solution = self.solution.sol
        subset = int(self.solution.nsteps * ratio)
        #print('driving frequency', self.model.driving_freq)
        period = 2 * np.pi /self.model.driving_freq#round(2 * np.pi /self.model.driving_freq,4)
        timestep = self.solution.tmax/int(self.solution.nsteps)
        #print('period = ',period)
        #print('timestep = ', timestep)
        dt = int(period/timestep)+1
        #print('dt = ', dt)
        fig = plt.figure()
        #https://stackoverflow.com/questions/14827650/pyplot-scatter-plot-marker-size
        plt.scatter(solution.y[0, subset::dt], solution.y[1, subset::dt], s=(pixel*72./fig.dpi)**2, marker="s")
        plt.xlabel('current', fontsize=18)
        plt.ylabel('d/dt current', fontsize=18)
        if name != 0:
        	plt.savefig(name, bbox_inches = "tight")
        fig.tight_layout()
        return fig
    
    def plot_phaseplane(self, ratio=0, pixel=1):
        solution = self.solution.sol
        subset = int(self.solution.nsteps * ratio)
        fig = plt.figure()
        plt.scatter(solution.y[0, subset:], solution.y[1, subset:], s=(pixel*72./fig.dpi)**2, marker="s")
        plt.xlabel('current', fontsize=18)
        plt.ylabel('d/dt current', fontsize=18)
        fig.tight_layout()
        return fig
        
    def plot_histogram(self):
        solution = self.solution.sol
        fig = plt.figure()
        hist, bin_edges = np.histogram(solution.y[0,:], density=True, bins=100)
        [x, y] = [bin_edges[:-1], hist]
        plt.plot(x, y) #'-ok', color='black'
        plt.plot([1.1*min(x), 1.1*max(x)], [0,0], color='black')
        plt.xlim((1.1*min(x), 1.1*max(x)))
        plt.xlabel('current', fontsize=18)
        plt.ylabel('pdf', fontsize=18)
        fig.tight_layout()
        return fig 
    
    def plot_potential(self,  xrange=np.linspace(-2,2,500)):
        x = xrange
        V = lambda x: -x*x/4 + x*x*x*x/8
        if self.model.potential_params != 0:
            #dV2 = lambda x:  0.5 * (-x + x*x*x) * (1 + a * np.cos(f*x))
            [a, f, s] = self.model.potential_params
            if s == 0:
                V2 = lambda x: (4* a*(-6+f*f*(-1+3*x*x))*np.cos(f*x) + 
                    f*x*(f**3*x*(-2+x*x) + 
                        4*a*(-6+f*f*(-1+x*x))*np.sin(f*x)))/(4*f**4)
            else:
                print('Plot unsymmetrical potential with only 2 minima')
                V2 = lambda x: s*x -x*x/4 + x*x*x*x/8
        
        fig = plt.figure()
        plt.plot(x, V(x))
        if self.model.potential_params != 0:
            plt.plot(x, V2(x))
        plt.xlabel('current', fontsize=18)
        plt.ylabel('potential', fontsize=18)
        #plt.savefig(model['name']+'poincare.png')
        fig.tight_layout()
        return fig
