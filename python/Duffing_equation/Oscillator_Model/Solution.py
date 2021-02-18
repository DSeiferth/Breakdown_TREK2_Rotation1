#
# Solution class
#
import numpy as np
import scipy.integrate
import os
import time

class Solution:
    """A  solution for osciallator equation for ion channels
    Parameters
    ----------
    model: using model class
        with model( gamma, initial_Condition)
    tmax: float
        integrates until it reaches tmax
        default value is 1
    nsteps: int
        number of integration steps
        default value is 1000
    """
    def __init__(self, model, tmax=1, nsteps=1000):
        self.model = model
        self.t_eval = np.linspace(0, tmax, nsteps)
        self.tmax = tmax
        self.nsteps = nsteps
        if model.driving_freq == 0:
            self.rand = np.random.normal(0,model.driving_amplitude,nsteps+1)

        self.sol = self.solver()

    def rhs(self, t, y):
        '''
        x'' = -gamma - dV(x)/dx + driving(t)
        y1 = x
        y2 = x'
        Right hand side of osciallator equation
        1st term : damping
        2nd term : derivative of potential
        3rd term : driving
        Parameters
        ----------
        t: time
        y: state vector [x, x']
        '''
        state = y
        gamma = self.model.gamma  # damping constant
        dt = self.tmax/self.nsteps
        
        dy_dt = [0, 0]  # [dqc, dq_p1, dq_p2]
        dy_dt[0] = state[1]
        if self.model.driving_freq > 0:
            # periodic driving
            driving = self.model.drivingterm(t)
        else:
            # stochastic driving
            step = int(t/dt)
            driving = self.rand[step]
        dy_dt[1] = -gamma * state[1] - self.model.dV_potential( state[0] ) + driving
        return dy_dt
    
    def solver(self):
        '''
        Runge-Kutta solver
        '''
        print('solving equation ...')
        start = time.time()
        step_func = self.rhs
        # intial condition
        self.y0 = self.model.IC
        self.sol = np.zeros(self.model.dim)

        sol = scipy.integrate.solve_ivp(
            fun=lambda t, y: step_func(t, y),
            t_span=[self.t_eval[0], self.t_eval[-1]],
            y0=self.y0,
            t_eval=self.t_eval
            #max_step=self.tmax / self.nsteps
        )
        self.sol = sol
        end = time.time()
        print('elapsed time for solving = ', end - start)
        return sol
    
    def output(self, folder):
        if not os.path.exists(folder):
            os.makedirs(folder)
        
        output = folder + '/time.txt'
        #np.savetxt(output, self.sol.t, delimiter=',')
        output = folder + '/current.txt'
        np.savetxt(output, self.sol.y[0])
        
        output = folder + '/parameters.txt'
        f = open(output, 'w')
        f.write('gamma \t {:f} \n'.format(self.model.gamma))
        f.write('driving_amplitude \t {:f} \n'.format(self.model.driving_amplitude))
        f.write('driving_freq \t {:f} \n'.format(self.model.driving_freq))
        f.write('initial_Cond[0] \t {:f} \n'.format(self.model.IC[0]))
        f.write('initial_Cond[1] \t {:f} \n'.format(self.model.IC[1]))
        if self.model.potential_params != 0:
            f.write('potential_params[0] \t {:f} \n'.format(self.model.potential_params[0]))
            f.write('potential_params[1] \t {:f} \n'.format(self.model.potential_params[1]))
        f.write('tmax {:d} \n'.format(self.tmax))
        f.write('nsteps {:d} \n'.format(self.nsteps))
        f.close()
