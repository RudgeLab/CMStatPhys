import numpy as np

def divide_RGB_log(parent,d1,d2,sigma):
    d1.log_red_protein = parent.log_red_protein + np.random.normal(0.0,sigma)
    d1.log_green_protein = parent.log_green_protein + np.random.normal(0.0,sigma)
    d1.log_blue_protein = parent.log_blue_protein + np.random.normal(0.0,sigma)
    
    d2.log_red_protein = parent.log_red_protein + np.random.normal(0.0,sigma)
    d2.log_green_protein = parent.log_green_protein + np.random.normal(0.0,sigma)
    d2.log_blue_protein = parent.log_blue_protein + np.random.normal(0.0,sigma)
    
    return parent,d1,d2

def divide_RGB_poisson(parent,d1,d2,lambd):
    d1.red_poisson = parent.red_poisson + np.random.poisson(lambd)
    d1.green_poisson = parent.green_poisson + np.random.poisson(lambd)
    d1.blue_poisson = parent.blue_poisson + np.random.poisson(lambd)
    
    d2.red_poisson = parent.red_poisson + np.random.poisson(lambd)
    d2.green_poisson = parent.green_poisson + np.random.poisson(lambd)
    d2.blue_poisson = parent.blue_poisson + np.random.poisson(lambd)
    
    return parent,d1,d2

def divide_binomial(parent,d1,d2):
    
    protein_division = np.random.binomial(parent.red_protein,0.5)
    
    d1.red_protein = protein_division
    d2.red_protein = (parent.red_protein-protein_division)
    
    return parent,d1,d2