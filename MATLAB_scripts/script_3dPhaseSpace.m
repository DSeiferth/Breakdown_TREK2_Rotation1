clc; clear all; close all; format compact;
addpath('../../0_Matlab_Functions/')
mV = -100;
downsampling = 1;
title_name = append('With noise ', num2str(mV),'mV')
title_name2 = '4-state Markov model (star)'

folder = append(num2str(mV), 'mV/');
%fname = append(num2str(mV), 'mV/', num2str(mV), '_mV_noNoise.txt');
fname = append(num2str(mV), 'mV/', num2str(mV), '_mV_Noise.txt');


table = readtable(fname);
conv =  table2array(table);

subset = 1e5
tau = 5
ndim = 3;
N_scatter = min(1e4, subset)-(ndim-1)*tau
subset = min(subset, size(conv,1))
x_t = conv(1:downsampling:subset*downsampling, 1);
disp('size of x_t')
disp(size(x_t))
mini = min(x_t)
maxi = max(x_t)
range = abs(maxi - mini)

plotname = append(title_name, '_tau_', num2str(tau),'.eps' )

Plot3dPhaseSpace(x_t,ndim, tau, subset, N_scatter, plotname, title_name, title_name2)    

