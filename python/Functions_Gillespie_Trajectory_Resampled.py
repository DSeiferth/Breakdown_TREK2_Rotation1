import numpy as np
from matplotlib import pyplot as plt
import math

def sum_up(matrix, row):
    res=0
    for i in range(len(matrix)):
        res=res+matrix[row][i]
        
    return res

def compare(vector, number):
    res=0
    for x in vector:
        if number<=vector[res]:
            return res
        else:
            res=res+1
    return -1

def rate_interval(matrix, row):
    vec=[]
    for col in range(1,len(matrix)+1):
        sum=0
        for i in range(col):
            sum=sum+matrix[row][i]
        res=sum/sum_up(matrix, row)
        vec.append(res)
    #print (vec)
    return vec


#initialisierung 
def init(matrix):
    sumUp=[]
    rateInterval=[]
    for row in range(len(matrix)):
        sumUp.append(sum_up(matrix,row))
        rateInterval.append(rate_interval(matrix,row))
    return [sumUp, rateInterval]
    #print(sumUp)
    #print(rateInterval)



def timestep(oldstate, t1, sumUp, rateInterval):
    dt=-np.log(np.random.random_sample() )/sumUp[oldstate]
    #print(dt)
    t1=t1+dt
    z1=np.random.random_sample()
    #print(z1)
    newstate=compare(rateInterval[oldstate], z1)
    #print(newstate)
    return [newstate, t1]

def simulation_traj(matrix, startState, T, sampling_step):
    [sumUp, rateInterval]=init(matrix)
    t=0.0
    state=startState
    #memory for trajectory
    TIME=[]
    states=[]

    while t<T:
        [state, t]=timestep(state, t, sumUp, rateInterval)
        TIME.append(t)
        states.append(state)
    prob=np.zeros(len(matrix))
    
    N = int(T/sampling_step)
    sample_states = -np.ones(N)
    sample_step = 0
    
    count_too_small = 0
    for i in range(len(TIME)-1):
        dt_gillespie = TIME[i+1]-TIME[i]
        n = int(dt_gillespie/sampling_step)
        if n < 1:
            count_too_small = count_too_small+1
        else:
            sample_states[sample_step:sample_step+n] = states[i]
            sample_step = sample_step +n
        prob[ states[i] ]=prob[ states[i] ] + dt_gillespie
    prob=1.0/T*prob
    print('number of Gillespie steps', len(states))
    print('number of sample steps ', sample_step)
    print('number of transitions not seen by sample step', count_too_small)
    return [prob, sample_states]
