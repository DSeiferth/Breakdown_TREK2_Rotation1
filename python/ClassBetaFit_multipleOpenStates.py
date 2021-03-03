import numpy as np
from matplotlib import pyplot as plt
import math

import time
import random
import networkx as nx
import os

from scipy.stats import norm
from scipy import signal
import Functions_Gillespie_Trajectory_Resampled as Gillespie
import Steady_State_Calculation_Spanning_Trees as auto

class ModelBetaFit:
    def __init__(self, T=1, samplingStep=2e-7, order=4, cutoff_fs=5e3, ResamplingFreq=200e3, expData=None, binning=100):
        self.__T = T
        self.__samplingStep = samplingStep
        self.__order = order
        self.__cutoff_fs = cutoff_fs
        self.__ResamplingFreq = ResamplingFreq
        ReSamplingStep = 1/ResamplingFreq
        self.__downsampling = int(ReSamplingStep/samplingStep)
        print('downsampling step = ', self.__downsampling )
        self.__exp_edges = binning
        if expData:
            self.__mV = expData[0]
            expFolder = expData[1]
            print(expFolder)
            fname_N = expFolder + str(int(cutoff_fs/1e3)) + 'kHz_' + str(self.__mV) + 'mV_N.txt'
            print('experimental data from: ', fname_N)
            fname_edges = expFolder +  str(int(cutoff_fs/1e3)) + 'kHz_' + str(self.__mV) + 'mV_edges.txt'
            self.__exp_N = np.loadtxt(fname_N);
            self.__exp_edges = np.loadtxt(fname_edges)
        
    def FIT_param(self, A, states, currentLevels, baselineNoise):
        self.__A = A
        self.__closed = states[0]
        self.__open = states[1]
        self.__closedLevel = currentLevels[0]
        self.__openLevel = currentLevels[1]
        self.__addNoise = baselineNoise[0]
        self.__convolutionNoise = baselineNoise[1]
        
    def GenerateGillespieTrajectory(self):
        G=auto.Matrix2Graph(self.__A)
        p=auto.steady_state(G)

        prob, sample_states, NoGillespie, NoMissed = Gillespie.simulation_traj(matrix=self.__A, startState=0, T=self.__T, sampling_step=self.__samplingStep)
        self.__NoGillespie = NoGillespie
        self.__NoMissed = NoMissed	
	
        sample_indistinguishable = sample_states.copy()
        # set open state to open currentLeve
        index = sample_states == -1
        sample_indistinguishable[index] = self.__openLevel[0]
        for i in range(len(self.__openLevel)):
            index = sample_states == self.__open[i]
            sample_indistinguishable[index] = self.__openLevel[i]

        # make state 0 and 1, 2, 3 indistinguishable
        for i in self.__closed:
            index = sample_states == i
            sample_indistinguishable[index] = self.__closedLevel[0]
        print('(analytical) steady state p = ',p)
        print('(simulated) Gillespie prob = ', prob)
        N = len(sample_states)
        print('number of sample states N = ', N)
        print('sample states min = ', sample_states.min())
        print('sample states max = ', sample_states.max())
        print('indistinguishable states min = ', sample_indistinguishable.min())
        print('indistinguishable states max = ', sample_indistinguishable.max())
        state0 = np.count_nonzero(sample_states == 0)
        state1 = np.count_nonzero(sample_states == 1)
        state2 = np.count_nonzero(sample_states == 2)
        state3 = np.count_nonzero(sample_states == 3)
        state4 = np.count_nonzero(sample_states == 4)
        error = np.count_nonzero(sample_states == -1)
        print('number of -1 states (error)', error, 'in percent', error/N)
        p_sample = [state0/N, state1/N, state2/N ,state3/N, state4/N ]
        print('p_sample = ', p_sample)
        self.__sample_indistinguishable = sample_indistinguishable
        print('unique values in sample_indistinguishable', np.unique(sample_indistinguishable))
        
        # Adding noise
        noise1 = np.random.normal(0, self.__addNoise,N)
        print('std. dev. of noise ', np.std(noise1))
        self.__signal_noise = sample_indistinguishable + noise1
        print()
        
    def Filter(self, plot=0, name=None, subset_ratio = 0.01):
        sampling_fs = 1/self.__samplingStep
        #print('sampling freq = ', sampling_fs)
        #print('sampling step = ', self.__samplingStep)

        # For analog filters, Wn is an angular frequency (e.g., rad/s).
        # For digital filters, Wn are in the same units as fs.
        #Type of output: numerator/denominator (‘ba’), pole-zero (‘zpk’), or second-order sections (‘sos’). 
        self.__sos = signal.bessel(self.__order, Wn=self.__cutoff_fs, btype='low', analog=False, output='sos', norm='phase', fs=sampling_fs)
        self.__y_sos = signal.sosfilt(self.__sos, self.__sample_indistinguishable) # no noise
        self.__y_sos_noise = signal.sosfilt(self.__sos, self.__signal_noise) # added noise
        
        if plot:
            subset = int(subset_ratio * self.__T/self.__samplingStep )
            f_size =  14
            fig, ax = plt.subplots()
            fig.suptitle('cut-off frequency ' + str(int(self.__cutoff_fs/1e3)) + 'kHz', fontsize=f_size) #str(mV)+'mV; ' +
            t = np.arange(0,self.__T, self.__samplingStep)
            ax.plot(t[0:subset], self.__sample_indistinguishable[0:subset], color='yellow', label='not filtered')
            ax.plot(t[0:subset], self.__y_sos[0:subset], '-', color='black', label='filtered')
            ax.set_xlabel('Time [s]', fontsize=f_size)
            ax.set_ylabel('Current [pA]', fontsize=f_size)
            ax.tick_params(axis='both', labelsize=f_size)
            ax.legend(prop={'size': 14})
            #ax.set_ylim([-1.1, 0.3])
            fig.tight_layout()
            if name:
                fig.savefig( name + str(int(self.__cutoff_fs/1e3)) + 'kHz_'+'timeseries_5states.pdf')
            plt.show()
    
    def NoiseLevelBeforeFilter(self):
        N = len(self.__sample_indistinguishable)
        print('length of sampled states = ', N)
        noise = np.random.normal(0, self.__addNoise,N)
        print('std. dev. of additive noise before filtering ', np.std(noise))
        # filter
        noise_sos = signal.sosfilt(self.__sos, noise)
        print('std. dev. of noise after filtering', np.std(noise_sos))
        #downsampling
        noise_sos = signal.sosfilt(self.__sos, noise)
        print('std. dev. of noise after filtering and downsampling', np.std(noise_sos[0::self.__downsampling]))
        print()
        
    def GenerateHistogram(self):
        binning = self.__exp_edges
        # Generate amplitude histograms
        [self.__hist_noNoiseFilter, self.__binEdges_noNoiseFilter]=np.histogram(self.__y_sos[0::self.__downsampling], density=True, bins=binning) # filtered
        [self.__hist_NoiseFilter, self.__binEdges_NoiseFilter]=np.histogram(self.__y_sos_noise[0::self.__downsampling], density=True, bins=binning) # filtered
        
        # Add noise with convolution
        rv1 = norm(loc = 0., scale = self.__convolutionNoise)
        dx1 = self.__binEdges_noNoiseFilter[2]-self.__binEdges_noNoiseFilter[1]
        x = np.arange(-10, 10, dx1)
        noise = rv1.pdf(x) 
        self.__hist_NoiseConvolution = signal.convolve(self.__hist_noNoiseFilter, noise, mode='same') / sum(noise)
    
    def PlotAmplitudeDistribution(self, log=1, limits = [[-75, 5], [1e-7, 1e0]],legend=0, exp=0, name=None, f_size=14):
        fig, ax = plt.subplots()
        #fig.suptitle('Gillespie trajectory sampled with %d kHz \n 4-pole Bessel filter, cut-off freq. %f kHz\n Re-sampling freq. %d kHz' %(int(sampling_fs/1e3) , int(ReSamplingFreq/1e3)), fontsize=f_size)
        ax.plot(self.__binEdges_noNoiseFilter[:-1], self.__hist_NoiseFilter, '--',  label='noise (added pre filter): $\sigma$ = %.2f' %self.__addNoise, linewidth=3)
        ax.plot(self.__binEdges_noNoiseFilter[:-1], self.__hist_NoiseConvolution, ':',  label='noise (convolution after filter): $\sigma$ = %.2f' %self.__convolutionNoise, linewidth = 3)

        if exp:
            ax.plot(self.__exp_edges, self.__exp_N, label = '%d mV recording' %self.__mV)
        ax.set_xlabel('Current [pA]', fontsize=f_size)
        ax.set_ylabel('PDF', fontsize=f_size)
        ax.set_xlim(limits[0])
        ax.tick_params(axis='both', labelsize=f_size)
        if log:
            ax.set_ylim(limits[1])
            ax.set_yscale('log')
            ax.plot([self.__openLevel, self.__openLevel], limits[1],label='current level of open state', color='black')
        if legend:
            ax.legend(prop={'size': 11}, loc='upper left')
        plt.tight_layout()
        if name:
            mV = self.__mV
            plt.savefig(name + str(int(self.__cutoff_fs/1e3)) + 'kHz_' +str(mV)+'_mV.pdf', bbox_inches='tight')
            if log:
                plt.savefig(name + str(int(self.__cutoff_fs/1e3)) + 'kHz_' +str(mV)+'_mV_log.pdf', bbox_inches='tight')
        plt.show()

    def WriteOutTrajectories(self, outFolder):
        downsampling = self.__downsampling #+str(self.__mV)+'_mV
        np.savetxt(outFolder + str(int(self.__cutoff_fs/1e3)) + 'kHz_noNoise.txt', self.__y_sos[0::downsampling])
        np.savetxt(outFolder + str(int(self.__cutoff_fs/1e3)) + 'kHz_Noise.txt', self.__y_sos_noise[0::downsampling])
        N_down = len(self.__sample_indistinguishable[0::downsampling] )
        traj =  self.__sample_indistinguishable[0::downsampling] + np.random.normal(0, self.__convolutionNoise,N_down)
        np.savetxt(outFolder + str(int(self.__cutoff_fs/1e3)) + 'kHz_NoFilter_BaselineNoise.txt', traj)
    
    @property
    def NoGillespie(self):
        print('Number of Gillespie events: ', self.__NoGillespie)
        return self.__NoGillespie
    @property
    def NoMissed(self):
        print('Number of Missed events after digitalizing Gillespie: ', self.__NoMissed)
        return self.__NoMissed
    @property
    def downsampling(self):
        return self.__downsampling
    @property
    def filtered_traj_noNoise(self):
        return self.__y_sos[0::self.__downsampling]
    @property
    def filtered_traj(self):
        return self.__y_sos_noise[0::self.__downsampling]
        
    @property
    def T(self):
        """
        length of simulated trajectory in [s]
        default = 1s
        """
        return self.__T
    @property
    def samplingStep(self):
        """
        sampling Step of Gillespie trajectory in [s] 
        used to discretize continuous traejctory
        default = 2e-7
        """
        return self.__samplingStep
    @property
    def order(self):
        """
        order of Bessel filter
        default order = 4
        """
        return self.__order
    @property
    def cutoff_fs(self):
        """
        cutoff frequency of Bessel filter
        default = 5e3 Hz
        """
        return self.__cutoff_fs
    @property
    def ResamplingFreq(self):
        """
        Resampling Frequency = sampling frequency of experimental set-up
        default = 200e3 Hz
        """
        return self.__ResamplingFreq

