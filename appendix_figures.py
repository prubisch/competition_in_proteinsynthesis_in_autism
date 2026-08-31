#standard python imports
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.style as style
import os

#import from .py files
from plotting_helper import *
from parameters import *
from model_variations import unlimited_r, limited_r, competitive_r

#custom colourbling friendly colourmaps
from cmcrameri import cm
palette = cm.batlow(np.linspace(0,0.9,6)) 

#check that directories exists for figures and simulation results
figures_path = 'figures/'
sim_path = 'simulation_results/'

models = [unlimited_r, limited_r, competitive_r]
R_nums = [Rs_unlim, Rs_lim, Rs_comp]


def run_sim_random_populations(model, r_range, n_runs): 


    results_population = np.zeros([n_runs, r_range.shape[0],3,5])
    
    for i in range(n_runs): 
        l_samples_rand = rng.gamma(mean_l, std_l, N)
        l_samples_rand = np.sort(l_samples_rand)
        
        b1 = np.nonzero(l_samples_rand < 50)[0][-1]
        b2 = np.nonzero(l_samples_rand < 100)[0][-1]
        b3 = np.nonzero(l_samples_rand < 150)[0][-1]
        b4 = np.nonzero(l_samples_rand < 200)[0][-1]

        s_1k = np.arange(0,b1+1,1, dtype = int)
        s_2k = np.arange(b1+1,b2+1,1, dtype = int)
        s_3k = np.arange(b2+1,b3+1,1, dtype = int)
        s_4k = np.arange(b3+1,b4+1,1, dtype = int)
        b_4k = np.arange(b4+1,N,1, dtype = int)


        #prepare parameter dictionary
        parameters = dict()

        parameters['α'] = α
        parameters['β'] = β/l_samples_rand

        parameters['γ'] = γ
        parameters['l'] = l_samples_rand
        parameters['lp'] = stiffness_fun(parameters['l'],stiffness_fac)
        parameters['Y_max'] = N
        parameters['N'] = N

        x0 = np.ones([N])*x_init


        for r in range(r_range.shape[0]):

            parameters['R'] = r_range[r]
            init_state =np.array([x0, np.zeros([N]),np.zeros([N]),np.zeros([N]),np.zeros([N]),np.zeros([N]),np.zeros([N])+parameters['R'],np.zeros([N])+parameters['Y_max']]) 

            res = simulation_wrapper(init_state,parameters,model, sim_time, dt, [True, sampling])
            for g_ind, g in enumerate([s_1k, s_2k, s_3k, s_4k, b_4k]): 
                results_population[i,r,0,g_ind]= np.mean(res[-1,2,g])
                results_population[i,r,1,g_ind]= np.std(res[-1,2,g])
                results_population[i,r,2,g_ind]= np.median(res[-1,2,g])
    
    
    return results_population


def appendix_figure_random_populations_vis_R_dep_y_prod(results, Rs_h): 
        
    layout = '''
        AAAAAA
        AAAAAA
        BBBBBB
        '''

    fig = plt.figure(figsize=(16, 6))

    height_ratios = [1,-0.2,1]
    width_ratios = [1,0.2,1,0.2,1,0.2]
    


    specs, gs = panel_specs(layout, fig=fig)
    gs.set_height_ratios(height_ratios)
    gs.set_width_ratios(width_ratios)
    gs.hspace=0.8
    axes = {}
    for letter in 'A':
        axes[letter] = ax = fig.add_subplot(specs[letter])
    #label_panels(axes, letters='A', postfix='', offset_left=1.2)
    for a in list('A'):
        axes[a].set_axis_off()

    mod_state = ['Unlimited R','Spatially restricted R', 'Competetive R']
    labels = ['<1k', '1-2k', '2-3k', '3-4k', '>4k']
    offset = 0.25

    for i in range(3): 
        ax = fig.add_subplot(gs[0,i*2])
        ax_st = fig.add_subplot(gs[2,i*2])
        ax.set_title(mod_state[i])
        
        
        for g in range(len(labels)): 
            #add example traces? 
            ax.plot(Rs_h[i], np.mean(results[i][:,:,0,g],axis = 0), label = labels[g], color = palette[g])
            ax.plot(Rs_h[i], np.mean(results[i][:,:,2,g],axis = 0), linestyle = ':')
            ax_st.plot(Rs_h[i], np.mean(results[i][:,:,1,g],axis = 0), label = labels[g], color = palette[g])
            
        if i ==2:
            ax.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left',  frameon = False)

        ax.set_xlabel('R')
        ax.set_ylabel('p')
        ax_st.set_xlabel('R')
        ax_st.set_ylabel('$\sigma(p)$')
        clear_axes(ax)
        clear_axes(ax_st)
        
    
    plt.savefig(figures_path + 'appendix_random_population_protein_yield_all_ribosome_conditions.pdf')

if __name__ == "__main__":
    
    if not os.path.isdir(figures_path):
        os.makedirs(figures_path)
    if not os.path.isdir(sim_path):
        os.makedirs(sim_path)
        
    #standard parameters
    file_ext = ''
    num_runs = 100
        
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
        print('Provide number of random population samples (100 is equivalent to provided data):')
        valid_num = False
        while not valid_num:
            ans_num = input()
            if ans_num.isdigit(): 
                num_runs = int(input)
                valid_num = True
                print('Warning long compute time needed! Are you sure you want to proceed? (y/n)')
                ans2 = input()
                if ans2 is in ['y','yes']: 
                    print('Running simulations with '+str(num_runs)+' random population samples.')
                else:
                    valid_num = False
                    print('Provide an alternative sampling size <'+ans_num)
            else: 
                print('Invalid input! Provide a numerical value.')
                
        if num_runs != 100: 
            file_ext = '_'+str(num_runs)+'_random_samples'
        for m_ind, m in enumerate(models): 
            data = run_sim_random_populations(m, R_nums[m_ind], num_runs)
            np.save(sim_path+f_names[m_ind]+file_ext+'.npy', data)
        
    print('Plotting Figure for random population sample statistics. (Appendix)')
    data_coll = []
    for m_ind, m in enumerate(models): 
        data_coll.append(np.load(sim_path+f_names[m_ind]+file_ext+'.npy'))
            
    appendix_figure_random_populations_vis_R_dep_y_prod(data_coll, R_nums)
        
        

