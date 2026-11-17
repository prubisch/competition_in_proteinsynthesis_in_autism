#standard python imports
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.style as style
import os

#import for data analytics of GO and CCDS data
from ftplib import FTP
import gzip
from scipy import stats

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
data_path = 'CCDS_data/'

models = [unlimited_r, limited_r, competitive_r]
R_nums = [Rs_unlim, Rs_lim, Rs_comp]
f_names = ['unlimited_R_effective_rates', 'limited_R_effective_rates', 'competitive_R_effective_rates']


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

    mod_state = ['Unlimited R','Spatially restricted R', 'Competitive R']
    labels = ['<1k', '1-2k', '2-3k', '3-4k', '>4k']
    offset = 0.25

    for i in range(3): 
        ax = fig.add_subplot(gs[0,i*2])
        ax_st = fig.add_subplot(gs[2,i*2])
        ax.set_title(mod_state[i])
        
        
        for g in range(len(labels)): 
            #add example traces? 
            ax.plot(Rs_h[i], np.mean(results[i][:,:,0,g],axis = 0), label = labels[g], color = palette[g])
            ax.plot(Rs_h[i], np.mean(results[i][:,:,2,g],axis = 0), linestyle = ':', color = palette[g])
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
    
    
######## functions for bio-informatics analytics of GO and CCDS data 
def retrieve_raw_data(species, to_retrieve, saved_f): 
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd('pub/CCDS/current_'+species+'/')
    with open(data_path+saved_f+'_raw_'+species+'.txt', 'wb') as fp:
        ftp.retrbinary('RETR '+to_retrieve, fp.write)
    fp.close()
    ftp.quit()
    
def retrieve_aa_data(species, to_retrieve, saved_f): 
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd('pub/CCDS/current_'+species+'/')
    
    #get the binary zip file
    with open('temp_file.faa.gz', 'wb') as f: 
        ftp.retrbinary('RETR '+to_retrieve, f.write)
    ftp.quit()
    f.close()
    
    #open the local zip file
    with gzip.open('temp_file.faa.gz', 'rb') as f, open(data_path+saved_f+'_raw_'+species+'.txt', 'wb') as fp: 
        fp.write(f.read())
    
    fp.close()
    f.close()
    #delete the local file
    os.remove('temp_file.faa.gz')
    
    
def retrieve_CCDS_ids_and_lengths_for_GO_genes(species, GO_gene_excel): 
    pframe = pd.read_csv(data_path+'CCDS_id_mapping_raw_'+species+'.txt', sep = '	')
    gframe = pd.read_excel(data_path+GO_gene_excel+'.xlsx')
    aframe = pd.read_excel(data_path+'CCDS_aa_sequence_'+species+'.xlsx')
    
    id_list = []
    length_list = []
    av_length_list = []
    sd_length_list = []
    
    #loop through the (SYN)GO list and retrieve the gene names
    for g_name in gframe['hgnc_symbol']: 
        #get the CCDS_ids from the CCDS_id-gene mapping file
        inds = np.where(pframe['gene'].str.lower()==g_name.casefold())[0]
        #check if there is an element in the index list, multiple will be stored in same column
        if inds.shape[0] == 0: 
            id_list.append([])
            length_list.append([np.nan])
            av_length_list.append(np.array([np.nan]))
            sd_length_list.append(np.array([np.nan]))
        else: 
            n_list = []
            l_list = np.array([])
            for i in inds: 
                #get the exact CCDS id for each match to the gene
                ccds_id = pframe['ccds_id'][i]
                n_list.append(ccds_id)
                #look up the length for the specific CCDS_id
                l_ind = np.where(aframe['CCDS_ID'] == ccds_id)[0]
                #check that ccds_id has a length
                if l_ind.shape[0] == 0: 
                    l_list = np.append(l_list, np.nan)
                else: 
                    l_list = np.append(l_list, aframe['#aa'][l_ind[0]])
                
            id_list.append(n_list)
            length_list.append(l_list)
            #calculate mean, check that the slice is not empty
            if np.sum(np.nan_to_num(l_list)) == 0.0: 
                av_length_list.append([np.nan])
                sd_length_list.append([np.nan])
            else: 
                av_length_list.append(np.nanmean(l_list))
                sd_length_list.append(np.nanstd(l_list))
    
    gframe.insert(len(gframe.keys()), 'CCDS_ID', id_list)
    gframe.insert(len(gframe.keys()), '#aa', length_list)
    gframe.insert(len(gframe.keys()), 'μ(#aa)', av_length_list)
    gframe.insert(len(gframe.keys()), 'σ(#aa)', sd_length_list)
    gframe.to_excel(data_path+GO_gene_excel+'_'+species+'.xlsx')
    
    
def extract_protein_length_matches(species, file_name): 
    csv_content = []
    
    aa_file = open(data_path+file_name+'_raw_'+species+'.txt')
    #per ccds_id
    lines = aa_file.read().split('>')
    for li in lines[1:]: 
        #extract ccds_id
        ccds_id, hs_num, aa_seq_p = li.split('|')
        #combined the full aa_seq into one string
        aa_seq = ''.join(aa_seq_p.split('\n')[1:])
        #get length
        csv_content.append([ccds_id,aa_seq,len(aa_seq)])
    
    
    df = pd.DataFrame(csv_content, columns=['CCDS_ID', 'AA_sequence', '#aa'])
    
    df.to_excel(data_path+file_name+'_'+species+'.xlsx', index = False)
    
    #close and delete the unnecessary raw file
    aa_file.close()
    os.remove(data_path+file_name+'_raw_'+species+'.txt')
    
    
def get_raw_data_files(species): 
    retrieve_raw_data(species,'CCDS.current.txt', 'CCDS_id_mapping')
    retrieve_aa_data(species, 'CCDS_protein.current.faa.gz', 'CCDS_aa_sequence')
    extract_protein_length_matches(species, 'CCDS_aa_sequence')

def check_for_overlaps(reference_frame, comparison_frame):
    
    reference_gene_list = reference_frame['hgnc_symbol'].to_list()
    compare_gene_list = comparison_frame['hgnc_symbol'].to_list()
    
    overlap_list = [gene for gene in reference_gene_list if gene in compare_gene_list]
    
    filtered_reference_frame = reference_frame[~reference_frame['hgnc_symbol'].isin(overlap_list)]
    filtered_comparison_frame = comparison_frame[~comparison_frame['hgnc_symbol'].isin(overlap_list)]
            
    return filtered_reference_frame, filtered_comparison_frame, overlap_list


def convert_list_of_strings_to_array(clist): 
    
    carray = np.array([])
    
    for item in clist: 
        carray = np.append(carray, np.fromstring(item.strip('[]'), dtype = float, sep = ' '))
    return carray


def plot_histo(GO_gene_excel1, GO_gene_excel2, label1, label2, remove_duplicates = True): 
    gframe = pd.read_excel(data_path+GO_gene_excel1+'.xlsx')
    rframe = pd.read_excel(data_path+GO_gene_excel2+'.xlsx')
    
    if remove_duplicates: 
        gframe.drop_duplicates(subset = ['hgnc_symbol'])
        rframe.drop_duplicates(subset = ['hgnc_symbol'])
    
    #ignoring overlapping genes
    l1 = convert_list_of_strings_to_array(gframe['#aa'])
    mean_l1 = gframe['μ(#aa)'].to_list()

    
    l2 =convert_list_of_strings_to_array(rframe['#aa'])
    mean_l2 = rframe['μ(#aa)'].to_list()
    
    #perform wilcoxon ranksum test
    statistics, pvalue = stats.ranksums(l1, l2, nan_policy = 'omit')
    
    
    filtered_gframe, filtered_rframe, overlap = check_for_overlaps(gframe, rframe)
    filtered_l1 = convert_list_of_strings_to_array(filtered_gframe['#aa'])
    filtered_l2 = convert_list_of_strings_to_array(filtered_rframe['#aa'])
    
    #perform wilcoxon ranksum test
    statistics_filtered, pvalue_filtered = stats.ranksums(filtered_l1, filtered_l2, nan_policy = 'omit')
    
    
    
    layout = '''
    AABB
    AABB
    '''

    fig = plt.figure(figsize=(12, 6))

    height_ratios = [1,0.1]
    width_ratios = [1,0.1,1,0.1]


    specs, gs = panel_specs(layout, fig=fig)
    gs.set_height_ratios(height_ratios)
    gs.set_width_ratios(width_ratios)
    gs.hspace=0.8
    axes = {}
    for letter in 'AB':
        axes[letter] = ax = fig.add_subplot(specs[letter])
    #label_panels(axes, letters='AB', postfix='', offset_left=1.2)
    for a in list('AB'):
        axes[a].set_axis_off()
        
    bin_width = 50
    binwin = np.arange(0,2000, bin_width)
        
    ax = fig.add_subplot(gs[0,0])
    ax.hist(l1, label = label1, color = palette[0], bins = binwin)
    ax.hist(l2, label = label2, color = palette[4], bins = binwin, alpha = 0.75)
    #ax.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left',  frameon = False)
    ax.set_xlabel('$l$ [#aa]')
    ax.set_ylabel('# genes')
    ax.set_title('All genes')
    ax.text(900,100, 'p = %.4e'%pvalue)
    clear_axes(ax)
    
    ax = fig.add_subplot(gs[0,2])
    ax.hist(filtered_l1, label = label1, color = palette[0], bins = binwin)
    ax.hist(filtered_l2, label = label2, color = palette[4], bins = binwin, alpha = 0.75)
    ax.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left',  frameon = False)
    ax.set_xlabel('$l$ [#aa]')
    ax.set_ylabel('# genes')
    ax.set_title('Unique genes')
    ax.text(900, 115, 'Excluded '+str(len(overlap))+' genes.')
    ax.text(900,100, 'p = %.4e'%pvalue_filtered)
    clear_axes(ax)
    
    
    plt.savefig(figures_path+'appendixfig_genelength_histograms'+label1+'_'+label2+'.pdf')
    

if __name__ == "__main__":
    
    if not os.path.isdir(figures_path):
        os.makedirs(figures_path)
    if not os.path.isdir(sim_path):
        os.makedirs(sim_path)
        
    #standard parameters
    file_ext = ''
    num_runs = 100
        
    print('Plot figures from provided simulation data? (y/n)')
    yn_valid = False
    while not yn_valid: 
        ans1 = input()
        if ans1 in ['y','yes','n','no']: 
            yn_valid = True
        else: 
            print('Invalid input! Valid input to indicate whether to plot figures from provided data: (y/n) or (yes/no)')
        
    if ans1 in ['y','yes']: 
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
                num_runs = int(ans_num)
                valid_num = True
                print('Warning long compute time needed! Are you sure you want to proceed? (y/n)')
                ans2 = input()
                if ans2 in ['y','yes']: 
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
            np.save(sim_path+f_names[m_ind]+'random_populations_'+file_ext+'.npy', data)
        
    print('Plotting Figure for random population sample statistics. (Appendix)')
    data_coll = []
    for m_ind, m in enumerate(models): 
        data_coll.append(np.load(sim_path+'random_populations_'+f_names[m_ind]+file_ext+'.npy'))
            
    appendix_figure_random_populations_vis_R_dep_y_prod(data_coll, R_nums)
    
    print('Download/Updating CCDS data? (y/n)')
    yn_valid = False
    while not yn_valid: 
        ans1 = input()
        if ans1 in ['y','yes','n','no']: 
            yn_valid = True
        else: 
            print('Invalid input! Valid input to indicate whether to download CCDS data: (y/n) or (yes/no)')
        
    if ans1 in ['y','yes']: 
        get_raw_data_files('human')
        get_raw_data_files('mouse')
        
    else: 
        if not os.path.isfile(data_path+'CCDS_aa_sequence_human.xlsx') and not os.path.isfile(data_path+'CCDS_aa_sequence_mouse.xlsx') and not os.path.isfile(data_path+'CCDS_id_mapping_raw_human.txt') and not os.path.isfile(data_path+'CCDS_id_mapping_raw_mouse.txt'): 
            print('CCDS data does not (fully) exist. Downloading data.')
            get_raw_data_files('human')
            get_raw_data_files('mouse')
    

    retrieve_CCDS_ids_and_lengths_for_GO_genes('human', 'syngo_1.3_process-in-the-postsynapse_annotations')
    retrieve_CCDS_ids_and_lengths_for_GO_genes('mouse', 'syngo_1.3_process-in-the-postsynapse_annotations')

    retrieve_CCDS_ids_and_lengths_for_GO_genes('human', 'syngo_1.3_metabolism_annotations')
    retrieve_CCDS_ids_and_lengths_for_GO_genes('mouse', 'syngo_1.3_metabolism_annotations')

    retrieve_CCDS_ids_and_lengths_for_GO_genes('human', 'syngo_1.3_synaptic-signaling_annotations')
    retrieve_CCDS_ids_and_lengths_for_GO_genes('mouse', 'syngo_1.3_synaptic-signaling_annotations')
    
    print('Plotting histograms and performing Wilcoxon rand sum test for SynGO sets.')
    plot_histo('syngo_1.3_synaptic-signaling_annotations_human', 'syngo_1.3_metabolism_annotations_human', 'synaptic signalling', 'metabolism')

    plot_histo('syngo_1.3_process-in-the-postsynapse_annotations_human', 'syngo_1.3_metabolism_annotations_human', 'postsynaptic processes', 'metabolism')
    
    

