import numpy as np
from matplotlib import pyplot as plt
import math
import time
import random
import networkx as nx
import os

#from scipy.stats import norm
#from scipy import signal
#import Functions_Gillespie_Trajectory_Resampled as Gillespie
#import Steady_State_Calculation_Spanning_Trees as auto
#import ClassBetaFit as BetaFit

def generateEventList(filtered_traj, threshold, step=5e-6, method='mean', closedLevel=0, openLevel=1):
    prev = filtered_traj[0]
    duration = []
    OpenClosed = []
    dur = 0
    idx_start = 0
    for i in range(1,len(filtered_traj)):
        if (prev > threshold and filtered_traj[i] > threshold) or (prev < threshold and filtered_traj[i] < threshold):
            dur = dur + step
        else:
            duration.append(dur)
            
            if method == 'mean':
                if idx_start == i :
                    meanCurrent = [i]
                else:
                    meanCurrent = np.mean(filtered_traj[idx_start:i])
                OpenClosed.append(meanCurrent)
            else:
                if prev > threshold: # closed
                    OpenClosed.append(closedLevel)
                elif prev < threshold: # open
                    OpenClosed.append(openLevel);
            prev = filtered_traj[i]
            idx_start = i
            dur = 0
    duration = np.array(duration)
    OpenClosed = np.array(OpenClosed)
    return duration, OpenClosed

def transformedDwellTime(x,x0):
    res = np.zeros(len(x))
    for i in range(len(x)):
        z = x[i]-x0
        res[i] = np.exp(z - np.exp(z))
    return res

def transformedDwelltimeStar(matrix, x):
    openTime = 1/(matrix[3][0] + matrix[3][1] + matrix[3][2] )
    slowTime = 1/matrix[0][3] 
    medTime = 1/matrix[1][3] 
    fastTime = 1/matrix[2][3] 
    
    slowArea = matrix[3][0] * openTime
    medArea = matrix[3][1] * openTime
    fastArea = matrix[3][2] * openTime
    
    slow = slowArea * transformedDwellTime(x,np.log(slowTime))
    med = medArea * transformedDwellTime(x,np.log(medTime))
    fast = fastArea * transformedDwellTime(x,np.log(fastTime))
    total = slow + med + fast
    return total, slow, med, fast
    
