import numpy as np

def DelayEmbedding(traj,ndim, tau, datcnt):
    delay = np.arange(0,(ndim-0)*tau, tau ) # [0 4 8]
    NoPhaseSpacePoints = datcnt-(ndim-1)*tau
    dataplot = np.zeros((NoPhaseSpacePoints, ndim))
    for i in range(NoPhaseSpacePoints):
        for j in range(ndim):
            d = delay[j]
            dataplot[i][j] = traj[i+d]
    return dataplot
