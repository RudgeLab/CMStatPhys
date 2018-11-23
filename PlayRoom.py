import LoadModifiedCellstates as LMC
import AddProteins as App
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
#ak

def load_data(nframes,lambd):
    
    cellstate = LMC.not_another_main(nframes,lambd)
    print "Loaded, adding radius"

    print "radius ready"
    #cellstate = LMC.load_unmod(1500)
    return cellstate

#--------------------------------------------------------------------
def calculate_sum_prot(cellstate,t):
    n = sum([cell.red_protein for id,cell in cellstate[t].iteritems()])
    return n

def get_R_max_t(cells_t):     #get max R

    R_max_t = 0.0
    for id,cell in cells_t.iteritems():
        if cell.r_dist > R_max_t:
            R_max_t = cell.r_dist
    return R_max_t

    
#def calculate_number_density(cellstates,N,t):
#    for id,cell in cellstates[t]:
#        sum_n_N = cell.redprot - N
def bin_check(cell_dist,r_bins,j):
    if cell_dist >= r_bins[j] and cell_dist < r_bins[j+1]:
        return True
    else:
        return False
    
#--------------------------------------------------------------------
#obtener curvas convergentes:
def obtain_convergent_curves(cellstate,t1,t2,nbins):
    n_norm = []
    for t in range(t1,t2):    
        
        cells_t = cellstate[t]
        
        R_max_t = get_R_max_t(cells_t)
                        
        N_prot = calculate_sum_prot(cellstate,t) #total num proteins in t
        
        n_prot = N_prot/(np.pi*(R_max_t**2))
        
        r_bins = np.linspace(0,1,nbins)
        
        n_i= np.array([])
        
        for j in range(0,len(r_bins)-1):
            
            bin_value = 0
            
            for id,cell in cells_t.iteritems():
                
                cell_dist = cell.r_dist/R_max_t
                
                if bin_check(cell_dist,r_bins,j) == True: 
                    
                    bin_value += cell.red_protein
                else:
                    bin_value += 0.0
                    
            area_i = np.pi*((r_bins[j+1])**2-(r_bins[j])**2)
    
            bin_value = bin_value/area_i
            
            n_i = np.append(n_i,bin_value)
            
        n_n = (n_i-n_prot)
        n_norm.append(n_n/max(n_i))
    return np.array(n_norm),r_bins
#--------------------------------------------------------------------
def calc_average_vol(cellstate):
    volume = []
    for t in range(len(cellstate)):
        for id,cell in cellstate[t].iteritems():
            r1 = cell.radius
            l1 = cell.length
            area1 = np.pi*r1**2 + 2*r1*l1
            volume.append(area1)
    volume = np.array(volume)
    avg_vol = np.average(volume)
    return avg_vol

def get_volume(cellstate_t):
    volume = []
    for id,cell in cellstate_t.iteritems():
        r1 = cell.radius
        l1 = cell.length
        area1 = np.pi*r1**2 + 2*r1*l1
        volume.append(area1)
    volume_t = np.sum(volume)
    return volume_t
            

def idea2(cellstate,t1,t2):    
    
    volume_t = np.array([get_volume(cellstate[t]) for t in range(t1,t2)])
    avg_vol = calc_average_vol(cellstate)

    print "average volume: ",avg_vol
    
    N_t = np.array([len(cellstate[t]) for t in range(t1,t2)])    
    
    R_max_t = np.array([get_R_max_t(cellstate[t]) for t in range(t1,t2)])
    
    R_t = np.sqrt(N_t*np.average(volume_t)/np.pi) 
    #R_t does not match R_max_t, missing space between cells
    
    
    return volume_t,avg_vol,N_t,R_max_t,R_t

def fix_Radii(volume_t,avg_vol,R_max_t,N_t):
    #volume expected by R_max_t:
    V_t = np.pi*(R_max_t**2)
    #we calculate the volume lost/cell
    Vv = (V_t-volume_t)/N_t #void volume per cell, converges to about 1
    #add new value to R_t as shown:
    R_t = np.sqrt(N_t*(avg_vol+Vv[-1])/np.pi)
    return R_t,Vv #correctamente se comapra con R_max_t


#--------------------------------------------------------------------
def obtain_Ik(N_t,t1,t2):
    dN = np.array([np.float(N_t[t+1] - N_t[t]) for t in range(t1,t2-2)])
    
    return dN
#--------------------------------------------------------------------
def obtain_n_t(cellstate,t1,t2):
    n_t = []
    for t in range(t1,t2):
        n_sum_t = calculate_sum_prot(cellstate,t)
        n_t.append(n_sum_t)
    return n_t
#--------------------------------------------------------------------
def add_grow_rate_i(cellstate,t):
    for id,cell in cellstate[t].iteritems():
        if t == 0:
            cell.mu = 0
        try:
            cellstate[t+1][id]
            r1 = cell.radius
            l1 = cell.length
            r2 = cellstate[t+1][id].radius
            l2 = cellstate[t+1][id].length
            area1 = np.pi*r1**2 + 2*r1*l1
            area2 = np.pi*r2**2 + 2*r2*l2
            
            mu_i = (1.0/area1)*(area2-area1)
            cell.mu = mu_i
        except:
            mu_i = cellstate[t-1][id].mu
            cell.mu = mu_i
    return cellstate[t]



def mu_bin_r(cellstate,t,nbins = 20):
    
    cells_t = cellstate[t]
    
    R_max_t = get_R_max_t(cells_t)
        
    r_bins = np.linspace(0,1,nbins)
    
    mu_r = np.array([])
    for j in range(0,len(r_bins)-1):
        bin_value = 0
        for id,cell in cells_t.iteritems():
            cell_dist = cell.r_dist/R_max_t
            if bin_check(cell_dist,r_bins,j) == True: 
                bin_value += cell.mu
            else:
                bin_value += 0.0
        mu_r = np.append(mu_r,bin_value)
        
    return mu_r,r_bins

def get_mu_t(cellstate,t1,t2):
    mu_t = np.array([get_mu_single_t(cellstate,t) for t in range(t1,t2)])
    return mu_t
        
    
def get_mu_single_t(cellstate,t):
    
    mu_total = 0
    cells_t = cellstate[t]
    
    for id,cell in cells_t.iteritems():
        mu_total += cell.mu
    mu_total = mu_total/len(cellstate[t])
    return mu_total
        
'''
def mu_bin_t(cellstate,t,nbins=20):
    
    cells_t = cellstate[t]
    
    T_max_t = len(cellstate)
        
    r_bins = np.linspace(0,1,nbins)
    
    mu_r = np.array([])
    for j in range(0,len(r_bins)-1):
        bin_value = 0
        for id,cell in cells_t.iteritems():
            cell_dist = cell.r_dist/R_max_t
            if bin_check(cell_dist,r_bins,j) == True: 
                bin_value += cell.mu
            else:
                bin_value += 0.0
        mu_r = np.append(mu_r,bin_value)
        
    return mu_r,r_bins
                    
'''
def plot_mu_r_t(r,t,mu_bins):
    ax = Axes3D(plt.gcf())
    ax.set_xlabel('r/R_max(t)')
    ax.set_ylabel('t/Nframes')
    ax.set_zlabel('mu/mu_max')
    ax.plot_wireframe(r,t,mu_bins)
    
def get_mu_r_t(cellstate,tbins = 897,nbins=50):
    for t in range(1,898):
        cellstate = add_grow_rate_i(cellstate,t)
    t_bins = np.linspace(0,1,tbins)
    r_bins = np.linspace(0,1,nbins)
    mu_bins = np.array([mu_bin_r(cellstate,t,nbins)[0] for t in range(1,tbins+1)])
    r_bins,t_bins = np.meshgrid(r_bins[:-1],t_bins)
    
    return r_bins,t_bins,mu_bins

def theoretical_N_t(cellstate,t1,t2):
    for t in range(t1,t2):
        cellstate[t] = add_grow_rate_i(cellstate,t)
    v_t,a_v,N_t,R_max_t,R_t = idea2(cellstate,t1,t2)
    R_t,Vv = fix_Radii(v_t,a_v,R_max_t,N_t)
    mu_t = get_mu_t(cellstate,t1,t2)
    int_mu = np.array([sum(mu_t[0:t]) for t in range(t1,t2)])
    N_tt = 1*np.exp(int_mu)
    R_tt = np.sqrt(v_t[0]/np.pi)*np.exp(0.5*int_mu)
    A_tt = v_t[0]*np.exp(int_mu)
    return N_tt,R_tt,A_tt,N_t,R_t,v_t,int_mu,Vv

cellstates = load_data(900,1)
N_tt,R_tt,A_tt,N_t,R_t,v_t,int_mu,Vv = theoretical_N_t(cellstates,0,899)
'''        
fig, ax = plt.subplots()
ax.plot(a, c, 'k--', label='Model length')
ax.plot(a, d, 'k:', label='Data length')
ax.plot(a, c + d, 'k', label='Total message length')

legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')
'''
#--------------------------------------------------------------------

#--------------------------------------------------------------------

#--------------------------------------------------------------------

#--------------------------------------------------------------------

#--------------------------------------------------------------------

#--------------------------------------------------------------------
    
#--------------------------------------------------------------------

#--------------------------------------------------------------------

#--------------------------------------------------------------------



    