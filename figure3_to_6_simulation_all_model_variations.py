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



def run_sim(model, r_range): 


    results_population = np.zeros([r_range.shape[0], int(sim_time/dt/sampling),8,N])

 
    #prepare parameter dictionary
    parameters = dict()
    
    parameters['α'] = α
    parameters['β'] = β/l_samples
    
    parameters['γ'] = γ
    parameters['l'] = l_samples
    parameters['lp'] = stiffness_fun(parameters['l'],stiffness_fac)
    parameters['Y_max'] = N
    parameters['N'] = N

    x0 = np.ones([N])*x_init
    

    for r in range(r_range.shape[0]):

        parameters['R'] = r_range[r]
        init_state =np.array([x0, np.zeros([N]),np.zeros([N]),np.zeros([N]),np.zeros([N]),np.zeros([N]),np.zeros([N])+parameters['R'],np.zeros([N])+parameters['Y_max']]) 

        res = simulation_wrapper(init_state,parameters,model, sim_time, dt, [True, sampling]) 
        results_population[r]= res
    
            
    return results_population

#some more analysis function which take the simulation results and calculate certain properties
#they are named according to what they calculate

def calc_means_over_t(results): 
    av1 = np.mean(results[:,s_1k], axis = 1)
    av2 = np.mean(results[:,s_2k], axis = 1)
    av3 = np.mean(results[:,s_3k], axis = 1)
    av4 = np.mean(results[:,s_4k], axis = 1)
    av5 = np.mean(results[:,b_4k], axis = 1)
    
    return [av1, av2,av3,av4,av5]
    
def calc_max_eff_rates_wr_R(results, rate_diff): 
    signed_maxes = np.zeros([results.shape[0],5])
    for r in range(results.shape[0]): 
        init_t = calc_means_over_t(results[r,:,rate_diff,:]) 
        elong_t = calc_means_over_t(results[r,:,rate_diff+1,:])
        for c in range(len(init_t)): 
            diff = init_t[c] - elong_t[c]
            ind_max = np.argmax(np.abs(diff))
            signed_maxes[r,c] = np.max(np.abs(diff)) * np.sign(diff[ind_max])
            
    return signed_maxes

def calc_tau_eff_rates_wr_R(results, rate_diff): 
    taus = np.zeros([results.shape[0],5])
    for r in range(results.shape[0]): 
        init_t = calc_means_over_t(results[r,:,rate_diff,:]) 
        elong_t = calc_means_over_t(results[r,:,rate_diff+1,:])
        for c in range(len(init_t)): 
            rate_t = np.abs(init_t[c]-elong_t[c])
            ind_max = np.argmax(rate_t)
            #this calculation assumes that we have an exponential decay, this is not always true
            #but close enough in approximation? 
            spacing = 1
            taus[r,c] = -dt*sampling*spacing/(np.log(rate_t[ind_max+spacing]/rate_t[ind_max]))
            #abs_max = np.amax(np.abs(rate_t))
            
            #half_point = np.nonzero(np.abs(rate_t[ind_max:]) < abs_max/2)[0][0]+ind_max 
            #taus[r,c] = half_point * dt*sampling* np.sign(rate_t[half_point]) #preserves relationship
            
    return taus
    
        
    

def calc_diff_wr_R(results,Rs_h, means = True):
    #just calculates the slope with regard to R
    #compute differences along correct axis
    #results need to have the shape [R, N]
    
    r_diff = np.diff(Rs_h)
    if means:
        mean_ys = np.mean(results, axis = 1)
        return np.diff(mean_ys)/r_diff
    
    else: 
        diff_r = np.diff(results, axis = 0)
        return diff_r/r_diff
    

def fig3_vis_R_dep_y_prod_ratio(results, Rs_h, indications_ind): 
        
    layout = '''
        AAAAAAAAAAA
        AAAAAAAAAAA
        BBBBBBBBBBB
        BBBBBBBBBBB
        CCCCCCCCCCC
        '''

    fig = plt.figure(figsize=(22, 10))

    height_ratios = [1,-0.2,1, -.2,1]
    width_ratios = [.5,0.2,.5,0.2,.5,0.2,.5,0.2,.5,.2,.5]
    


    specs, gs = panel_specs(layout, fig=fig)
    gs.set_height_ratios(height_ratios)
    gs.set_width_ratios(width_ratios)
    gs.hspace=0.8
    axes = {}
    for letter in 'ABC':
        axes[letter] = ax = fig.add_subplot(specs[letter])
    label_panels(axes, letters='ABC', postfix='', offset_left=1.2)
    for a in list('ABC'):
        axes[a].set_axis_off()

    mod_state = ['Unlimited R','Spatially restricted R', 'Competetive R']
    r_labels = ['I', 'II', 'III', 'IV']
    offset = 0.25

    for i in range(3): 
        ax = fig.add_subplot(gs[0,i*4:(i*4)+3])
        ax.set_title(mod_state[i])
        av1 = np.mean(results[i][:,-1,2,s_1k], axis = 1)
        av2 = np.mean(results[i][:,-1,2,s_2k], axis = 1)
        av3 = np.mean(results[i][:,-1,2,s_3k], axis = 1)
        av4 = np.mean(results[i][:,-1,2,s_4k], axis = 1)
        av5 = np.mean(results[i][:,-1,2,b_4k], axis = 1)
        
        

        ax.plot(Rs_h[i], av1, label = '<1k', color = palette[0])
        ax.plot(Rs_h[i], av2, label = '1-2k', color = palette[1])
        ax.plot(Rs_h[i], av3, label = '2-3k', color = palette[2])
        ax.plot(Rs_h[i], av4, label = '3-4k', color = palette[3])
        ax.plot(Rs_h[i], av5, label = '>4k', color = palette[4])
        miny, maxy = ax.get_ylim()
        diffmm = maxy - miny
        for k in range(len(indications_ind[i])):
            ax.plot([Rs_h[i][indications_ind[i][k]],Rs_h[i][indications_ind[i][k]]], [miny, maxy], linestyle = ':', color = 'grey')
            ax.text(Rs_h[i][indications_ind[i][k]],maxy-offset * diffmm,r_labels[k], horizontalalignment = 'center')
        if i ==2:
            ax.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left',  frameon = False)

        ax.set_xlabel('R')

        ax.set_ylabel('p')
        miny, maxy = ax.get_ylim()
        clear_axes(ax)
        
        
        ax = fig.add_subplot(gs[2,i*4:(i*4)+3])
        ax.plot(Rs_h[i][1:], calc_diff_wr_R(results[i][:,-1, 2, -1], Rs_h[i], means = False)/calc_diff_wr_R(results[i][:,-1,2,0], Rs_h[i], means = False),label = '$\\frac{max}{min}$', color = palette[0])
        ax.plot(Rs_h[i][1:], calc_diff_wr_R(results[i][:,-1,2,b_4k], Rs_h[i])/calc_diff_wr_R(results[i][:,-1,2,s_1k], Rs_h[i]), label = '$\\frac{>4k}{<1k}$', color = palette[2])
        ax.plot(Rs_h[i][1:], calc_diff_wr_R(results[i][:,-1,2,long_ind], Rs_h[i])/calc_diff_wr_R(results[i][:,-1,2,short_ind], Rs_h[i]), label = '$\\frac{>2k}{<2k}$', color = palette[4])
        miny, maxy = ax.get_ylim()
        diffmm = maxy - miny
        for k in range(len(indications_ind[i])):
            ax.plot([Rs_h[i][indications_ind[i][k]],Rs_h[i][indications_ind[i][k]]], [miny, maxy], linestyle = ':', color = 'grey')
            ax.text(Rs_h[i][indications_ind[i][k]],maxy-offset * diffmm,r_labels[k], horizontalalignment = 'center')
        

        ax.set_xlabel('R')

        ax.set_ylabel('$\\rho_{i,j}(R)$')
        ax.plot(Rs_h[i][1:], Rs_h[i][1:]* 0, linestyle = ':', color = 'grey')
        if i ==2:
            ax.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left',  frameon = False)
        clear_axes(ax)
        
        
        ax = fig.add_subplot(gs[4,i*4])
        ax.plot(Rs_h[i][1:], calc_diff_wr_R(results[i][:,-1, 2, -1], Rs_h[i], means = False),label = '$max$', color = palette[0])
        ax.plot(Rs_h[i][1:], calc_diff_wr_R(results[i][:,-1,2,b_4k], Rs_h[i]), label = '$>4k$', color = palette[2])
        ax.plot(Rs_h[i][1:], calc_diff_wr_R(results[i][:,-1,2,long_ind], Rs_h[i]), label = '$>2k$', color = palette[4])
        miny, maxy = ax.get_ylim()
        diffmm = maxy - miny
        for k in range(len(indications_ind[i])):
            ax.plot([Rs_h[i][indications_ind[i][k]],Rs_h[i][indications_ind[i][k]]], [miny, maxy], linestyle = ':', color = 'grey')
            ax.text(Rs_h[i][indications_ind[i][k]],maxy-offset * diffmm,r_labels[k], horizontalalignment = 'center')

            
        ax.set_xlabel('R')

        ax.set_ylabel('$\\frac{dp_{i}}{dR}$')
        ax.plot(Rs_h[i][1:], Rs_h[i][1:]* 0, linestyle = ':', color = 'grey')
        clear_axes(ax)
        
        ax = fig.add_subplot(gs[4,i*4+2])
        ax.plot(Rs_h[i][1:], calc_diff_wr_R(results[i][:,-1,2,0], Rs_h[i], means = False), label = '$i = max, j = min$', color = palette[0], linestyle = '-.')
        ax.plot(Rs_h[i][1:], calc_diff_wr_R(results[i][:,-1,2,s_1k], Rs_h[i]), label = '$i >4k, j <1k$', color = palette[2], linestyle = '-.')
        ax.plot(Rs_h[i][1:], calc_diff_wr_R(results[i][:,-1,2,short_ind], Rs_h[i]), label = '$i >2k, j <2k$', color = palette[4], linestyle = '-.')
        miny, maxy = ax.get_ylim()
        diffmm = maxy - miny
        for k in range(len(indications_ind[i])):
            ax.plot([Rs_h[i][indications_ind[i][k]],Rs_h[i][indications_ind[i][k]]], [miny, maxy], linestyle = ':', color = 'grey')
            ax.text(Rs_h[i][indications_ind[i][k]],maxy-offset * diffmm,r_labels[k], horizontalalignment = 'center')
        
        
        ax.set_xlabel('R')
        ax.set_ylabel('$\\frac{dp_{j}}{dR}$')
        ax.plot(Rs_h[i][1:], Rs_h[i][1:]* 0, linestyle = ':', color = 'grey')
        if i ==2:
            ax.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left',  frameon = False)
        clear_axes(ax)
    

    plt.savefig(figures_path + 'fig3_protein_yield_all_ribosome_conditions_ratios.pdf')
    
    
def fig4_vis_diff_effective_rates(results_all, ribo_ind, plot_till = int(dt*sampling*sim_time)): 
   
    layout = '''
        AAAAAA
        AAAAAA
        BBBBBB
        BBBBBB
        CCCCCC
        '''

    fig = plt.figure(figsize=(16, 8))

    height_ratios = [1,-0.2,1,-0.2,1]
    width_ratios = [1,0.2,1,0.2,1,0.2]


    specs, gs = panel_specs(layout, fig=fig)
    gs.set_height_ratios(height_ratios)
    gs.set_width_ratios(width_ratios)
    gs.hspace=0.8
    axes = {}
    for letter in 'ABC':
        axes[letter] = ax = fig.add_subplot(specs[letter])
    label_panels(axes, letters=['A: Unlimited R', 'B: Spatially restricted, limited R', 'C: Competitive R'], postfix='', offset_left=1.2)
    for a in list('ABC'):
        axes[a].set_axis_off()

    ribo_state = ['I','II', 'III']
    classes = ['<1k', '1-2k', '2-3k', '3-4k', '>4k']
    

    for j in range(3): 
        for i in range(3): 
            ax = fig.add_subplot(gs[j*2,i*2])
            if j == 0: 
                ax.set_title(ribo_state[i])
            #calculating the means of the rates we read out, 
            averages_alpha = calc_means_over_t(results_all[j][ribo_ind[j][i],:,3,:])
            averages_beta = calc_means_over_t(results_all[j][ribo_ind[j][i],:,4,:])
            for c in range(len(classes)): 
                ax.plot(np.arange(0,sim_time,dt*sampling)[:plot_till], averages_alpha[c][:plot_till]-averages_beta[c][:plot_till], label = classes[c], color = palette[c])
            ax.plot(np.arange(0,sim_time,dt*sampling)[:plot_till],np.arange(0,sim_time,dt*sampling)[:plot_till]*0,color = 'grey', linestyle = ':')
            ax.set_xlabel('t')
            ax.set_ylabel('$α_{eff,i}-β_{eff,i}$')
            if i == 2 and j == 0: 
                ax.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left',  frameon = False)
            clear_axes(ax)

    plt.savefig(figures_path+'fig4_diff_effective_rates_double_all_models.pdf')
    
    
    
def fig5_vis_resource_abundance(results_all, ribo_ind, plot_till = int(dt*sampling*sim_time)): 
    
    layout = '''
        AABBCC
        AABBCC
        '''

    fig = plt.figure(figsize=(16, 4))

    height_ratios = [1,-0.2]
    width_ratios = [1,0.4,1,0.4,1,0.2]


    specs, gs = panel_specs(layout, fig=fig)
    gs.set_height_ratios(height_ratios)
    gs.set_width_ratios(width_ratios)
    gs.hspace=0.8
    axes = {}
    for letter in 'ABC':
        axes[letter] = ax = fig.add_subplot(specs[letter])
    #label_panels(axes, letters='ABC', postfix='', offset_left=1.2)
    for a in list('ABC'):
        axes[a].set_axis_off()

    ribo_state = ['I','II', 'III']
    line_styles = ['solid','dotted', 'dashed']
    mod_state = ['Unlimited R','Spatially restricted, limited R', 'Competetive R']
    
    for j in range(3): 
        ax = fig.add_subplot(gs[0,j*2])
        ax2 = ax.twinx()
        for i in range(3): 
            ax.plot(np.arange(0,sim_time,dt*sampling)[:plot_till], results_all[j][ribo_ind[j][i],:plot_till,6,0], color = palette[0], label = ribo_state[i], linestyle = line_styles[i])
            ax2.plot(np.arange(0,sim_time,dt*sampling)[:plot_till], results_all[j][ribo_ind[j][i],:plot_till,7,0], color = palette[2], label = ribo_state[i], linestyle = line_styles[i])
        
        ax.plot(np.arange(0,sim_time,dt*sampling)[:plot_till],np.arange(0,sim_time,dt*sampling)[:plot_till]*0,color = 'grey', linestyle = ':')
        ax2.plot(np.arange(0,sim_time,dt*sampling)[:plot_till],np.arange(0,sim_time,dt*sampling)[:plot_till]*0,color = 'grey', linestyle = ':')
        ax.set_xlabel('t')
        ax2.set_xlabel('t')
        miny, maxy = ax.get_ylim()
        ax.set_ylim(0,maxy)
        miny, maxy = ax2.get_ylim()
        ax2.set_ylim(0,maxy)
        ax.set_title(mod_state[j])
        ax.set_ylabel('$r$', color = palette[0])
        ax2.set_ylabel('$y$', color = palette[2])
        ax.tick_params(axis='y', colors=palette[0])
        ax2.tick_params(axis='y', colors=palette[2])
        if i == 2 and j == 2: 
                ax.legend(bbox_to_anchor=(1.2, 1.0), loc='upper left',  frameon = False)
        clear_axes(ax)
        clear_second_axes(ax2)

    plt.savefig(figures_path+'fig5_resource_abundances.pdf')
    
    
def fig6_comp_y_prod_for_theory_visulation(results, Rs_h): 
    
    layout = '''
        AA
        AA
        '''

    fig = plt.figure(figsize=(8, 4))

    height_ratios = [1,-0.2]
    width_ratios = [1,0.1]
    


    specs, gs = panel_specs(layout, fig=fig)
    gs.set_height_ratios(height_ratios)
    gs.set_width_ratios(width_ratios)
    gs.hspace=0.8
    axes = {}
    for letter in 'A':
        axes[letter] = ax = fig.add_subplot(specs[letter])
    #label_panels(axes, letters='ABC', postfix='', offset_left=1.2)
    for a in list('A'):
        axes[a].set_axis_off()
 
    ax = fig.add_subplot(gs[0,0])
    av1 = np.mean(results[:,-1,2,s_1k], axis = 1)
    av2 = np.mean(results[:,-1,2,s_2k], axis = 1)
    av3 = np.mean(results[:,-1,2,s_3k], axis = 1)
    av4 = np.mean(results[:,-1,2,s_4k], axis = 1)
    av5 = np.mean(results[:,-1,2,b_4k], axis = 1)

    ax.plot(Rs_h, av1, label = '<1k', color = palette[0])
    ax.plot(Rs_h, av2, label = '1-2k', color = palette[1])
    ax.plot(Rs_h, av3, label = '2-3k', color = palette[2])
    ax.plot(Rs_h, av4, label = '3-4k', color = palette[3])
    ax.plot(Rs_h, av5, label = '>4k', color = palette[4])
    
    ax.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left',  frameon = False)
    ax.set_xlabel('R')
    ax.set_ylabel('p')
    clear_axes(ax)
    
    plt.savefig(figures_path + 'fig6_protein_yield_competitive_y_and_r.pdf')
    
    
    
#TODO: add a user input whether the data should be simulated or the figures simply replotted
    
#run the simulation and safe the data in .npy files to reload later for data visualisation
#this is just a safety net should the server time out during running the simulations or plotting


f_names = ['unlimited_R_effective_rates', 'limited_R_effective_rates', 'competitive_R_effective_rates']

 
    
#reload the saved data in case the server died and you dont want to rerun all simulations


#these are the ribosome indixes for which the effective rates and resource abundanced are computed
r_to_plot_ul = [1,10,40]
r_to_plot_l = [3,15,40]
r_to_plot_d = [8,23,40]




if __name__ == "__main__":
    
    if not os.path.isdir(figures_path):
        os.makedirs(figures_path)
    if not os.path.isdir(sim_path):
        os.makedirs(sim_path)
        
    print('Plot figures from provided input?')
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
        unlim_R_pop = run_sim(unlimited_r, Rs_unlim)
        np.save(sim_path+'unlimited_R_effective_rates.npy', unlim_R_pop)

        lim_R_pop = run_sim(limited_r, Rs_lim)
        np.save(sim_path+'limited_R_effective_rates.npy', lim_R_pop)

        comp_R_pop = run_sim(competitive_r, Rs_comp)
        np.save(sim_path+'competitive_R_effective_rates.npy', comp_R_pop)
        
    #plot data
    unlim_R_pop = np.load(sim_path+'unlimited_R_effective_rates.npy')
    lim_R_pop = np.load(sim_path+'limited_R_effective_rates.npy')
    comp_R_pop = np.load(sim_path+'competitive_R_effective_rates.npy')
    
    print('Figure 3')
    fig3_vis_R_dep_y_prod_ratio([unlim_R_pop,lim_R_pop, comp_R_pop], [Rs_unlim, Rs_lim, Rs_comp],[r_to_plot_ul, r_to_plot_l, r_to_plot_d])
    print('Figure 4')
    fig4_vis_diff_effective_rates([unlim_R_pop,lim_R_pop, comp_R_pop],[r_to_plot_ul, r_to_plot_l, r_to_plot_d], plot_till = 100)
    print('Figure 5')
    fig5_vis_resource_abundance([unlim_R_pop,lim_R_pop, comp_R_pop],[r_to_plot_ul, r_to_plot_l, r_to_plot_d], plot_till = 100)
    print('Base for figure 6 (Paper figures postprocessed with inkscape).')
    fig6_comp_y_prod_for_theory_visulation(comp_R_pop, Rs_comp, safe = True)
