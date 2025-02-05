import numpy as np


#rates
rng = np.random.RandomState(seed = 5)
γ = 0.05
β = 0.1
α = 0.1

#standard number of mRNAs, initial mRNA concentration per mRNA and ribosome range
N = 100
x_init = 1
Rs = np.linspace(0.001,200,50)

#mRNA length descriptors
mean_l, std_l = 2, 50


#sampling the lengths
l_samples = rng.gamma(mean_l, std_l, N)
l_samples = np.sort(l_samples)

#compute mRNA groups by size
b1 = np.nonzero(l_samples < 50)[0][-1]
b2 = np.nonzero(l_samples < 100)[0][-1]
b3 = np.nonzero(l_samples < 150)[0][-1]
b4 = np.nonzero(l_samples < 200)[0][-1]

s_1k = np.arange(0,b1+1,1, dtype = int)
s_2k = np.arange(b1+1,b2+1,1, dtype = int)
s_3k = np.arange(b2+1,b3+1,1, dtype = int)
s_4k = np.arange(b3+1,b4+1,1, dtype = int)
b_4k = np.arange(b4+1,N,1, dtype = int)

short_ind = np.concatenate((s_1k, s_2k), dtype = 'int')
long_ind = np.concatenate((s_3k, s_4k, b_4k), dtype = 'int')

#combinations from the reTASEP paper
l_combs = np.array([[100,500],[500,2500],[100,2500]])


#length dependent stiffness for reinitiation assumption
def R_fernandes(length, stiffness_mRNA, xr): 
    #might need to restrict to the 0 and 1 range... 
    return np.sqrt(2*length) * (xr + (1-xr)*stiffness_mRNA)**0.5

def stiffness_fun(length,k): 
    return length**(-k)

#factor how stiffness affect reinition rate
stiffness_fac = 0.5



#general simulation_wrapper
def simulation_wrapper(state0, params, update,  elong_fun,sim_time, dt, record = [True, 100]): 
    
    states = np.zeros([int(sim_time/dt/record[1])+1, state0.shape[0], state0.shape[1]])
    next_t = state0
    
    for t in range(int(sim_time/dt)):
        
        if record[0]:
            if t % record[1] == 0:
                states[int(t/record[1])+1] = next_t
        next_t = update(params, next_t, dt, elong_fun)
    states[-1] = next_t
        
    return states[1:]