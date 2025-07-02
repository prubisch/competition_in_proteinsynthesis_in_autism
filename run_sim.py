import os
import sys
import numpy as np
from parameters import l_samples,s_1k,s_2k,s_3k,s_4k,b_4k


bngl_dir = 'bngl_files/'
vizdata_dir = 'viz_data/seed_001/'
observable_dir = 'simulation_results/'

every_n = 1e4 #how often we produce viz data

ITERATIONS = 1e7

dt = 1e-6


def run_simulation(cond): 
    
    MCELL_PATH = os.environ.get('MCELL_PATH', '')
    sys.path.append(os.path.join(MCELL_PATH, 'lib'))

    import mcell as m


    viz_output = m.VizOutput(
        output_files_prefix = vizdata_dir+cond+'/Scene',
        every_n_timesteps= every_n
        )


    model = m.Model()

    model.add_viz_output(viz_output)

    model.load_bngl(
        bngl_dir+cond+'.bngl', 
        observables_path_or_file=observable_dir+cond+'.gdat')

    param_dict = m.bngl_utils.load_bngl_parameters(bngl_dir+cond+'.bngl')

    for count in model.counts:
        count.every_n_timesteps = every_n

    model.config.use_bng_units = False

    model.config.total_iterations = ITERATIONS
    model.config.time_step = dt


    model.config.partition_dimension = 1

    model.initialize()

    model.export_viz_data_model()

    model.run_iterations(ITERATIONS)

    model.end_simulation()
    
def write_sample_bngl(l,N): 
    #define the different l distributions here
    file_path = bngl_dir+'template_'+str(l)+'_'+str(N)+'.bngl'
    if os.path.exists(file_path):
        return file_path
    else:
        exp_file = open(bngl_dir+'example.bngl', 'r')
        new_file = open(file_path, 'w')
        if l[0] == '2': 
            l_dist = np.array([l[1],l[2]], dtype = int)
        if l[0] == 'dist': 
            l_dist = l_samples
            
        
        for line in exp_file: 
            if '%%%%insert numbers' in line: 
                for l_i in l_dist:
                    new_file.write('    nucl_'+str(l_i)+' '+str(N)+'\n')
                
            elif '%%%%insert nucl def' in line: 
                mid = ''
                for l_i in l_dist:
                    mid = mid+'~'+str(l_i)
                new_file.write('    nucl(l,r,R,i~0~1,t~0~1,s'+str(mid)+')\n')
                    
            elif '%%%%insert protein def' in line: 
                mid = ''
                for l_i in l_dist:
                    mid = mid+'~'+str(l_i)
                new_file.write('    P(s'+mid+')\n')
                
            elif '%%%%insert P release' in line: 
                for l_i in l_dist:
                    new_file.write('    @CP:P(s~'+str(l_i)+') P_i\n')
                
            
            elif '%%%%insert mRNA def' in line:
                for l_i in l_dist:
                    mid = ''
                    for x in range(l_i-2): #because initiation and termination are 2 parts as well 
                        mid = mid+'.nucl(R,t~0,i~0,s~'+str(l_i)+',l!'+str(x)+',r!'+str(x+1)+')'
                    new_file.write('@CP:nucl(R,t~0,i~1,s~'+str(l_i)+',r!0,l)'+mid+'.nucl(R,t~1,i~0,s~'+str(l_i)+',l!'+str(x+1)+',r) nucl_'+str(l_i)+'\n')
                    
            elif '%%%%insert degration equation' in line:
                for l_i in l_dist:
                    mid = ' + Y' * (l_i-1) 
                    new_file.write('    P(s~'+str(l_i)+') -> Y'+mid+' deg_r'+'\n')
            
            elif '%%%%insert reaction equation' in line:
                for l_i in l_dist:
                    new_file.write('    nucl(R!1, t~1, s~'+str(l_i)+').Ribo(N!1) + Y -> P(s~'+str(l_i)+')+ nucl(R, t~1, s~'+str(l_i)+')+ Ribo(N) term_r \n')
                
            elif '%%%%insert observ' in line:
                for l_i in l_dist: 
                    new_file.write('    Molecules p_'+str(l_i)+' P(s~'+str(l_i)+')\n')
                for l_i in l_dist: 
                    new_file.write('    Molecules r_initiated_'+str(l_i)+' nucl(R!1, i~1, s~'+str(l_i)+').Ribo(N!1)\n')
                for l_i in l_dist: 
                    new_file.write('    Molecules r_elongating_'+str(l_i)+' nucl(R!1, t~0, i~0, s~'+str(l_i)+').Ribo(N!1)\n')
                for l_i in l_dist: 
                    new_file.write('    Molecules r_terminating_'+str(l_i)+' nucl(R!1, t~1, i~0, s~'+str(l_i)+').Ribo(N!1)\n')
                    
            else:
                new_file.write(line)
        exp_file.close()
        new_file.close()
    return file_path
    
    
    
    

def check_and_write_bngl(R, Y, l,N):
    #the different distributions of mRNAs are predefined in the example file 
    #doing this automatically would mean a lot more inserting than just exchanging numbers
    #we will abuse chatGPT as a finite state machine to get those without too much typing
    condition = str(R)+'ribos_'+str(Y)+'aminos_'+str(l)+'_'+str(N)
    file_path = bngl_dir+condition+'.bngl'
    if os.path.exists(file_path):
        return condition
    else:
        template_path = bngl_dir+'template_'+str(l)+'_'+str(N)+'.bngl'
        if not os.path.exists(template_path):
            write_sample_bngl(l,N)
        exp_file = open(template_path, 'r')
        new_file = open(file_path, 'w')
        
        for line in exp_file: 
            if not '%%%%insert' in line: 
                #need to check whether line break is already in it. 
                new_file.write(line)
            elif '%%%%insert ribo' in line: 
                new_file.write('    Ribo_i '+str(R)+'\n')
                
            elif '%%%%insert Y' in line: 
                new_file.write('    Y_i '+str(Y)+'\n')
        
        
        exp_file.close()
        new_file.close()
    return condition
#get input from the "comand line" to know which simulation to start

if __name__ == "__main__":
    #parameters provided by the .sh script are: 
    #l: length description we want to test
    #N : number of mRNAs
    #R: number of ribosomes
    #Y: number of aminoacids

    sim_type = sys.argv[1]
    n_mRNA = int(sys.argv[2])
    n_ribo = int(sys.argv[3])
    n_amino = int(sys.argv[4])
    #these are only needed when sim_type ==2
    l_specs1 = int(sys.argv[5])
    l_specs2 = int(sys.argv[6])
    if sim_type == '2': 
        species_lengths = ['2',l_specs1,l_specs2]
    elif sim_type =='dist':
        species_lengths = sim_type
    else:
        raise Exception('Only 2 species or distribution setting allowed.') 
                            
    condition = check_and_write_bngl(n_ribo, n_amino, species_lengths, n_mRNA)
    run_simulation(condition)