import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from model_variations import *
from parameters import *

#plot the parameter exploration results
from cmcrameri import cm
from plotting_helper import *
figures_path = 'figures/'
sim_path = 'simulation_results/'
c_map = cm.berlin

import os

#check that the directories exists for both the figures and the simulation results


model_list = [ulR_ulY, ulR_lY, ulR_cY, lR_ulY, lR_lY, lR_cY, cR_ulY, cR_lY, cR_cY]


def slope_calc(α_h, β_h , model): 
    
    R1, R2 = [50,100]

 
    #prepare parameter dictionary
    parameters = dict()
    
    parameters['α'] = α_h
    parameters['β'] = β_h/l_samples
    
    parameters['γ'] = γ #γ_h
    parameters['l'] = l_samples
    parameters['lp'] = stiffness_fun(parameters['l'],stiffness_fac)
    parameters['Y_max'] = N
    parameters['N'] = N


    x0 = np.ones([N])*x_init
    init_state =np.array([x0, np.zeros([N]),np.zeros([N])]) 
    
    parameters['R'] = R1
    y1 = simulation_wrapper(init_state,parameters,model,R_fernandes, sim_time, dt, [True, sampling])[-1, 2,:]
    

    parameters['R'] = R2
    y2 = simulation_wrapper(init_state,parameters,model,R_fernandes, sim_time, dt, [True, sampling])[-1, 2,:]
    
    slope_here = (np.mean(y1[long_ind] - y2[long_ind])/ np.mean(y1[short_ind] - y2[short_ind]))
    
    return slope_here


def exploration_plot(filename, colorbar_minmax): 
    
    data = np.load(sim_path+filename+'.npy')
    
    data = data.reshape([alpha_vals.shape[0], beta_vals.shape[0]])
    
    layout = '''
        AA
        AA
        '''

    fig = plt.figure(figsize=(22, 10))

    height_ratios = [1,0.2]
    width_ratios = [1,0.2]
    


    specs, gs = panel_specs(layout, fig=fig)
    gs.set_height_ratios(height_ratios)
    gs.set_width_ratios(width_ratios)
    gs.hspace=0.8
    axes = {}
    for letter in 'A':
        axes[letter] = ax = fig.add_subplot(specs[letter])
    for a in list('A'):
        axes[a].set_axis_off()


    
    ax = fig.add_subplot(gs[0,0])
    pic = ax.imshow(data, interpolation = 'none', vmin = colorbar_minmax[0], vmax = colorbar_minmax[1], cmap = c_map)
    

    cbar = fig.colorbar(pic, ax=ax, label = 'ρ')
    ax.set_xticks(np.arange(beta_vals.shape[0])[::5])
    ax.set_xticklabels(np.round(beta_vals,3)[::5])
    ax.tick_params("x", rotation=45)
    ax.set_xlabel('β')
    ax.set_yticks(np.arange(alpha_vals.shape[0])[::5])
    ax.set_yticklabels(np.round(alpha_vals,3)[::5])
    ax.set_ylabel('α')

    
         
    plt.savefig(figures_path + 'fig2_small_fig'+filename+'.pdf')


def run_parallel(model, filename): 
    res = Parallel(n_jobs = -1)(delayed(slope_calc) (alpha, beta, model) for alpha in alpha_vals for beta in beta_vals)
    np.save(sim_path+filename+'.npy',res)


lims = [-0.5,0.5]

f_names = ['unlimitedR_unlimitedY_rho', 'unlimitedR_limitedY_rho','unlimitedR_competitiveY_rho', 'limitedR_unlimitedY_rho','limitedR_limitedY_rho','limitedR_competitiveY_rho','competitiveR_unlimitedY_rho','competitiveR_limitedY_rho','competitiveR_competitiveY_rho']

    
    
    
if __name__ == "__main__":
    
    if not os.path.isdir(figures_path):
        os.makedirs(figures_path)
    if not os.path.isdir(sim_path):
        os.makedirs(sim_path)
        
    #standard parameters
    file_ext = ''
    num_param = 50
    alpha_vals = np.linspace(0.001, 0.5,num_param)
    beta_vals = np.linspace(0.001, 0.5,num_param)
        
    print('Plot figures from provided data? (y/n)')
    yn_valid = False
    while not yn_valid: 
        ans1 = input()
        if ans1 is in ['y','yes','n','no']: 
            yn_valid = True
        else: 
            print('Invalid input! Valid input to indicate whether to plot figures from provided data: (y/n) or (yes/no)')
        
    if ans1 is in ['y','yes']: 
        #check if preload data exists
        rerun = False
        for f in f_names: 
            if not os.path.isfile(sim_path+f+'.npy'): 
                print('Provided data does not exists.')
                rerun = True
    else: 
        rerun = True

    if rerun: 
        print('Provide number of parameter values for the simulation (50 is equivalent to provided data):')
        valid_num = False
        while not valid_num:
            ans_num = input()
            if ans_num.isdigit(): 
                num_param = int(input)
                valid_num = True
                print('Warning long compute time needed! Are you sure you want to proceed? (y/n)')
                ans2 = input()
                if ans2 is in ['y','yes']: 
                    print('Running simulations with '+str(num_param)+' values for α and β.')
                else:
                    valid_num = False
                    print('Provide an alternative sampling size <'+ans_num)
            else: 
                print('Invalid input! Provide a numerical value.')
         
        #run the simulations
        if num_param != 50: 
            file_ext = '_'+str(num_param)+'_values_'
        alpha_vals = np.linspace(0.001, 0.5,num_param)
        beta_vals = np.linspace(0.001, 0.5,num_param)
        for m_ind, m in enumerate(model_list): 
            run_parallel(m, f_names[m_ind]+file_ext)
            
    print('Plotting single subfigures.')
    for f in f_names: 
        exploration_plot(f+file_ext, lims)
    print('Figure 2 of the paper is constructed from the single subfigures with inkscape.')
