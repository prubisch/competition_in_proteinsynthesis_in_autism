import numpy as np

model_descriptions = ['unlimited_R_unlimited_Y','unlimited_R_limited_Y','unlimited_R_competitive_Y','limited_R_unlimited_Y','limited_R_limited_Y','limited_R_competitive_Y','competitive_R_unlimited_Y','competitive_R_limited_Y','competitive_R_competitive_Y',]

def ulR_ulY(parameters, state_var, dt):

    
    x_vec = state_var[0]   #free mRNA concentrations
    xr_vec = state_var[1]  # bound mRNA ribosome concentrations
    y_vec = state_var[2]   #protein concentrations
    r = parameters['R']  #unlimited, constant
    y_free = parameters['Y_max'] # unlimited, constant


    new_x_vec = x_vec + dt * (-parameters['α']  * r * x_vec  + parameters['β']* xr_vec *y_free)
    new_xr_vec = xr_vec + dt * (parameters['α'] * r *  x_vec - parameters['β'] * xr_vec*y_free)
    new_y_vec = y_vec + dt * (parameters['β']  * xr_vec*y_free - parameters['γ'] * y_vec)


    return new_x_vec, new_xr_vec, new_y_vec


def ulR_lY(parameters, state_var, dt):

    
    x_vec = state_var[0]   #free mRNA concentrations
    xr_vec = state_var[1]  # bound mRNA ribosome concentrations
    y_vec = state_var[2]   #protein concentrations
    r = parameters['R']  #unlimited, constant
    y_free = parameters['Y_max']/(parameters['N']/2) - y_vec #limited, non competitive


    new_x_vec = x_vec + dt * (-parameters['α']  * r * x_vec  + parameters['β']* xr_vec *y_free)
    new_xr_vec = xr_vec + dt * (parameters['α'] * r *  x_vec - parameters['β'] * xr_vec*y_free)
    new_y_vec = y_vec + dt * (parameters['β']  * xr_vec*y_free - parameters['γ'] * y_vec)


    return new_x_vec, new_xr_vec, new_y_vec


def ulR_cY(parameters, state_var, dt):

    
    x_vec = state_var[0]   #free mRNA concentrations
    xr_vec = state_var[1]  # bound mRNA ribosome concentrations
    y_vec = state_var[2]   #protein concentrations
    r = parameters['R']  #unlimited, constant
    y_free = parameters['Y_max'] - np.sum(y_vec) # limited, competitve


    new_x_vec = x_vec + dt * (-parameters['α']  * r * x_vec  + parameters['β']* xr_vec *y_free)
    new_xr_vec = xr_vec + dt * (parameters['α'] * r *  x_vec - parameters['β'] * xr_vec*y_free)
    new_y_vec = y_vec + dt * (parameters['β']  * xr_vec*y_free - parameters['γ'] * y_vec)


    return new_x_vec, new_xr_vec, new_y_vec



def lR_ulY(parameters, state_var, dt):

    
    x_vec = state_var[0]   #free mRNA concentrations
    xr_vec = state_var[1]  # bound mRNA ribosome concentrations
    y_vec = state_var[2]   #protein concentrations
    r = parameters['R']/(parameters['N']/2) - xr_vec  #limited, uncompetitive
    y_free = parameters['Y_max']  #unlimited, constant


    new_x_vec = x_vec + dt * (-parameters['α']  * r * x_vec  + parameters['β']* xr_vec *y_free)
    new_xr_vec = xr_vec + dt * (parameters['α'] * r *  x_vec - parameters['β'] * xr_vec*y_free)
    new_y_vec = y_vec + dt * (parameters['β']  * xr_vec*y_free - parameters['γ'] * y_vec)


    return new_x_vec, new_xr_vec, new_y_vec

def lR_lY(parameters, state_var, dt):

    
    x_vec = state_var[0]   #free mRNA concentrations
    xr_vec = state_var[1]  # bound mRNA ribosome concentrations
    y_vec = state_var[2]   #protein concentrations
    r = parameters['R']/(parameters['N']/2) - xr_vec   #limited, uncompetitive
    y_free = parameters['Y_max']/(parameters['N']/2) - y_vec  #limited, uncompetitive


    new_x_vec = x_vec + dt * (-parameters['α']  * r * x_vec  + parameters['β']* xr_vec *y_free)
    new_xr_vec = xr_vec + dt * (parameters['α'] * r *  x_vec - parameters['β'] * xr_vec*y_free)
    new_y_vec = y_vec + dt * (parameters['β']  * xr_vec*y_free - parameters['γ'] * y_vec)


    return new_x_vec, new_xr_vec, new_y_vec

def lR_cY(parameters, state_var, dt):

    
    x_vec = state_var[0]   #free mRNA concentrations
    xr_vec = state_var[1]  # bound mRNA ribosome concentrations
    y_vec = state_var[2]   #protein concentrations
    r = parameters['R']/(parameters['N']/2) - xr_vec  #limited, uncompetitive
    y_free = parameters['Y_max'] - np.sum(y_vec) #limited, competitive


    new_x_vec = x_vec + dt * (-parameters['α']  * r * x_vec  + parameters['β']* xr_vec *y_free)
    new_xr_vec = xr_vec + dt * (parameters['α'] * r *  x_vec - parameters['β'] * xr_vec*y_free)
    new_y_vec = y_vec + dt * (parameters['β']  * xr_vec*y_free - parameters['γ'] * y_vec)


    return new_x_vec, new_xr_vec, new_y_vec


def cR_ulY(parameters, state_var, dt):

    
    x_vec = state_var[0]   #free mRNA concentrations
    xr_vec = state_var[1]  # bound mRNA ribosome concentrations
    y_vec = state_var[2]   #protein concentrations
    r = parameters['R'] - np.sum(xr_vec)  #limited, competitive
    y_free = parameters['Y_max'] #unlimited, noncompetitive


    new_x_vec = x_vec + dt * (-parameters['α']  * r * x_vec  + parameters['β']* xr_vec *y_free)
    new_xr_vec = xr_vec + dt * (parameters['α'] * r *  x_vec - parameters['β'] * xr_vec*y_free)
    new_y_vec = y_vec + dt * (parameters['β']  * xr_vec*y_free - parameters['γ'] * y_vec)


    return new_x_vec, new_xr_vec, new_y_vec

def cR_lY(parameters, state_var, dt):

    
    x_vec = state_var[0]   #free mRNA concentrations
    xr_vec = state_var[1]  # bound mRNA ribosome concentrations
    y_vec = state_var[2]   #protein concentrations
    r = parameters['R'] - np.sum(xr_vec)  #limited, competitive
    y_free = parameters['Y_max']/(parameters['N']/2) - y_vec #limited, noncompetitive


    new_x_vec = x_vec + dt * (-parameters['α']  * r * x_vec  + parameters['β']* xr_vec *y_free)
    new_xr_vec = xr_vec + dt * (parameters['α'] * r *  x_vec - parameters['β'] * xr_vec*y_free)
    new_y_vec = y_vec + dt * (parameters['β']  * xr_vec*y_free - parameters['γ'] * y_vec)


    return new_x_vec, new_xr_vec, new_y_vec

def cR_cY(parameters, state_var, dt):

    
    x_vec = state_var[0]   #free mRNA concentrations
    xr_vec = state_var[1]  # bound mRNA ribosome concentrations
    y_vec = state_var[2]   #protein concentrations
    r = parameters['R'] - np.sum(xr_vec)  #limited, competitive
    y_free = parameters['Y_max'] - np.sum(y_vec) #limited, competitive


    new_x_vec = x_vec + dt * (-parameters['α']  * r * x_vec  + parameters['β']* xr_vec *y_free)
    new_xr_vec = xr_vec + dt * (parameters['α'] * r *  x_vec - parameters['β'] * xr_vec*y_free)
    new_y_vec = y_vec + dt * (parameters['β']  * xr_vec*y_free - parameters['γ'] * y_vec)


    return new_x_vec, new_xr_vec, new_y_vec


#these are the model variants and there returns for gathering effective initiation and elongation rates
#these include all three models with the lcY (limited-competitive Y) condition
def unlimited_r(parameters, state_var, dt):

    
    x_vec = state_var[0]   #free mRNA concentrations
    xr_vec = state_var[1]  # bound mRNA ribosome concentrations
    y_vec = state_var[2]   #protein concentrations
    r = parameters['R']  #unlimited, constant
    y_free = parameters['Y_max'] - np.sum(y_vec) # limited, competitve
    
    effective_init = parameters['α']  * r * x_vec 
    effective_elong = parameters['β']* xr_vec *y_free
    effective_deg = parameters['γ'] * y_vec


    new_x_vec = x_vec + dt * (-effective_init + effective_elong)
    new_xr_vec = xr_vec + dt * (effective_init - effective_elong)
    new_y_vec = y_vec + dt * (effective_elong - effective_deg)


    return new_x_vec, new_xr_vec, new_y_vec, effective_init, effective_elong, effective_deg, r*np.ones(parameters['N']), y_free *np.ones(parameters['N'])

def limited_r(parameters, state_var, dt):

    
    x_vec = state_var[0]   #free mRNA concentrations
    xr_vec = state_var[1]  # bound mRNA ribosome concentrations
    y_vec = state_var[2]   #protein concentrations
    r = parameters['R']/(parameters['N']/2) - xr_vec  #limited, uncompetitive
    y_free = parameters['Y_max'] - np.sum(y_vec) #limited, competitive


    effective_init = parameters['α']  * r * x_vec 
    effective_elong = parameters['β']* xr_vec *y_free
    effective_deg = parameters['γ'] * y_vec


    new_x_vec = x_vec + dt * (-effective_init + effective_elong)
    new_xr_vec = xr_vec + dt * (effective_init - effective_elong)
    new_y_vec = y_vec + dt * (effective_elong - effective_deg)


    return new_x_vec, new_xr_vec, new_y_vec, effective_init, effective_elong, effective_deg, r, y_free *np.ones(parameters['N'])

def competitive_r(parameters, state_var, dt):

    
    x_vec = state_var[0]   #free mRNA concentrations
    xr_vec = state_var[1]  # bound mRNA ribosome concentrations
    y_vec = state_var[2]   #protein concentrations
    r = parameters['R'] - np.sum(xr_vec)  #limited, competitive
    y_free = parameters['Y_max'] - np.sum(y_vec) #limited, competitive


    effective_init = parameters['α']  * r * x_vec 
    effective_elong = parameters['β']* xr_vec *y_free
    effective_deg = parameters['γ'] * y_vec


    new_x_vec = x_vec + dt * (-effective_init + effective_elong)
    new_xr_vec = xr_vec + dt * (effective_init - effective_elong)
    new_y_vec = y_vec + dt * (effective_elong - effective_deg)


    return new_x_vec, new_xr_vec, new_y_vec, effective_init, effective_elong, effective_deg, r*np.ones(parameters['N']), y_free *np.ones(parameters['N'])
