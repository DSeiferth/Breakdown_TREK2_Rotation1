import numpy as np
from pyentrp import entropy as ent

def resampling(duration, OpenClosed, noiseLevel, samplingStep=5e-6, display=1):
    resampled = [];
    zero_count = 0;
    for idx in range(len(duration)):
        dur = duration[idx];
        no = np.floor(dur/samplingStep);
        if no == 0:
            zero_count = zero_count + 1
        else:
            sdata = np.ones(int(no)) * OpenClosed[idx];
            resampled =  np.concatenate((resampled, sdata))
    if display == 1:
        print('length of event List')
        print(len(duration))
        print('size of resampled vector')
        print(len(resampled))
        print('number of zero-count events')
        print(zero_count)
    # add noise
    if noiseLevel:
        resampled = resampled + np.random.normal(0, noiseLevel,len(resampled))
    return resampled


def surrogateAnalysisRandomShuffle(duration, OpenClosed, subset, noiseLevel, samplingStep, Delay=1, Nsurr=10):
    #print(len(duration))
    #print(len(OpenClosed))
    #print(noiseLevel)
    resampledEventList = resampling(duration[0:subset], OpenClosed[0:subset], noiseLevel, samplingStep)
    H_EvenList = ent.permutation_entropy(resampledEventList, order=5, delay=Delay, normalize=True)
    print('noise level ', noiseLevel)
    print('resampled event List H = ', H_EvenList)
    entropy_vector = np.zeros(Nsurr)
    NoEvents = len(duration)
    for i in range(Nsurr):
        permutation_vector = np.random.permutation(NoEvents)
        duration_perm = duration[permutation_vector]
        OpenClosed_perm = OpenClosed[permutation_vector]
        resampled_perm = resampling(duration_perm, OpenClosed_perm, noiseLevel, samplingStep, display=0)
        H_perm = ent.permutation_entropy(resampled_perm, order=5, delay=Delay, normalize=True)
        entropy_vector[i] = H_perm
    min_surr = min(entropy_vector)
    max_surr = max(entropy_vector)
    print('minimal entropy of surrogates: H = ', min_surr)
    print('maximal entropy of surrogates: H = ', max_surr)
    if H_EvenList < min_surr:
        print('Deterministic clasification')
    return H_EvenList, min_surr, max_surr
