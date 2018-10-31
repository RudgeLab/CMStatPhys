from RunAndCompress import GetSubDir
import DivisionDistribution as DD
import cPickle 
import os
import numpy as np

    
def load2pickles(fname,j): #loads pickle j and j+1
    data_1 = cPickle.load(open(fname%(j))) 
    data_2 = cPickle.load(open(fname%(j+1))) 

    cellstate_1 = data_1['cellStates'] 
    lineage_1 = data_1['lineage']
    cellstate_2 = data_2['cellStates'] 
    lineage_2 = data_2['lineage']
    
    return cellstate_1,lineage_1,cellstate_2,lineage_2

def loadPickle_lite(fname,j): #loads pickle j and j+1
    data_1 = cPickle.load(open(fname%(j))) 

    cellstate_1 = data_1['cellStates'] 
    lineage_1 = data_1['lineage']

    return cellstate_1,lineage_1

def add_proteins_2(cellstate_1,cellstate_2,lineage_2,t,sigma = 0, lambd = 0): #cellstate[t],cellstate[t+1],lineage[t+1]
    
    #------------- dividing module
    for id,parent in cellstate_1.iteritems():
        if t == 0:

            parent.red_protein = 0

        try:

            cellstate_2[id] #divides?
            cellstate_2[id].red_protein = parent.red_protein + np.random.poisson(lambd)

            
        except KeyError:

            dids = [key for key, value in lineage_2.iteritems() if value == id]
            did1 = dids[0]
            did2 = dids[1]

            #------Binomial division
            cellstate_1[id],cellstate_2[did1],cellstate_2[did2] = DD.divide_binomial(parent,cellstate_2[did1],cellstate_2[did2])
        #print '+++++++++++next cell'
    return cellstate_1,cellstate_2

def find_zero(cellstate):
    x_0,y_0,z_0 = 0.0,0.0,0.0
    for id,cell in cellstate.iteritems():
        x_0 += cell.pos[0]
        y_0 += cell.pos[1]
        z_0 += cell.pos[2]
    return x_0/len(cellstate),y_0/len(cellstate),z_0/len(cellstate)

def add_radius_angle_area(cellstate):
    
    x_0,y_0,z_0 = find_zero(cellstate)
    
    for id,cell in cellstate.iteritems():
        
        x_cell,y_cell,z_cell = cell.pos[0],cell.pos[1],cell.pos[2]
        cell.r_dist =  np.sqrt(((x_cell-x_0)**2)+((y_cell-y_0)**2)+((z_cell-z_0)**2))
        #angle
        angle_r = np.sqrt(x_cell**2 + y_cell**2 + z_cell**2)
        
        if y_cell-y_0 >= 0.0 and cell.r_dist != 0.0:

            cell.phi = np.arccos((x_cell-x_0)/cell.r_dist)
        elif y_cell-y_0 < 0.0 and cell.r_dist != 0.0:
            cell.phi = 2*np.pi - np.arccos((x_cell-x_0)/cell.r_dist)
        else:
            cell.phi = 0.0
        #area    
        r = cell.radius
        l = cell.length
        cell.area = np.pi*r**2 + 2*r*l
        
    return cellstate

def add_ndiv(cellstate_1,cellstate_2):
    try:
        cellstate_1[1].ndiv = 1
        cellstate_2[1].ndiv = 1
    except:
        a = 0
        
    for id,cell in cellstate_2.iteritems():
        try: #if doesn't divide register the new amount of cells in the current cellstate
            cellstate_1[id] 
            cell.ndiv = cellstate_1[id].ndiv
            #must be number of cells when the cell first appears
        except KeyError: #if it divides register the amount of cells at time of division
            cell.ndiv = len(cellstate_2)
        
    return cellstate_1,cellstate_2

def bin_check(cell_dist,r_bins,j):
    if cell_dist >= r_bins[j] and cell_dist < r_bins[j+1]:
        return True
    else:
        return False
    
def calculate_sum_prot(cellstate):
    n = sum([cell.red_protein for id,cell in cellstate.iteritems()])
    return n

def get_R_max_t(cells_t):     #get max R

    R_max_t = 0.0
    for id,cell in cells_t.iteritems():
        if cell.r_dist > R_max_t:
            R_max_t = cell.r_dist
    return R_max_t   

def get_max_area(cells_t):     #get max R

    max_area = 0.0
    for id,cell in cells_t.iteritems():
        if cell.area > max_area:
            max_area = cell.area
    return max_area   

def obtain_convergent_curves(cellstate,nbins=20):
    
    cells_t = cellstate
    
    R_max_t = get_R_max_t(cells_t)
                    
    N_prot = calculate_sum_prot(cells_t) #total num proteins in t
    try:
        n_prot = N_prot/(np.pi*(R_max_t**2))
    
    except ZeroDivisionError:
        cell = cellstate[1]
        r1 = cell.radius
        l1 = cell.length
        area1 = np.pi*r1**2 + 2*r1*l1

        n_prot = N_prot/area1
    
    r_bins = np.linspace(0,1,nbins)
    
    n_i= np.array([])
    
    for j in range(0,len(r_bins)-1):
        
        bin_value = 0
        
        for id,cell in cells_t.iteritems():
            
            if R_max_t != 0.0:
                cell_dist = cell.r_dist/R_max_t
            else:
                cell_dist = 1.0
            
            if bin_check(cell_dist,r_bins,j) == True: 
                
                bin_value += cell.red_protein
            else:
                bin_value += 0.0
                
        area_i = np.pi*((r_bins[j+1])**2-(r_bins[j])**2)

        bin_value = bin_value/area_i
        
        n_i = np.append(n_i,bin_value)
        
    return n_i,n_prot,r_bins

def obtain_convergent_curves_branch(branch,R_max_t,nbins=20):
    
    cells_t = branch
                        
    max_area = get_max_area(cells_t)
            
    r_bins = np.linspace(0,1,nbins)
    
    n_i= np.array([])
    
    for j in range(0,len(r_bins)-1):
        
        bin_value = 0
        area_bin = 0.0
        
        for id,cell in cells_t.iteritems():
            
            if R_max_t != 0.0:
                cell_dist = cell.r_dist/R_max_t
            else:
                cell_dist = 1.0
            
            if bin_check(cell_dist,r_bins,j) == True: 
                
                bin_value += cell.red_protein
                area_bin += cell.area/max_area
                
            else:
                
                bin_value += 0.0
                
        if area_bin != 0:
            bin_value = bin_value/area_bin
        
        n_i = np.append(n_i,bin_value)
        
    return n_i,r_bins
'''
datafolders = []
root = "/Users/Medina/cellmodeller"
rootdir = root+"/data"
    
startframe = 0      
    
datafolders,datafiles,folders = GetSubDir(rootdir)
i = 0

nframes = 700
lambd = 1.0
nbins = 25
repetitions = 3
path_to_write = root+"/convergent_curves"
if not os.path.isdir(path_to_write):
    os.makedirs(path_to_write)
    
for simulation in datafiles:
    
    for repet in range(1,repetitions+1):
        
        print "Repetition: ", repet
        print 'Loading and running '+ datafolders[i]
        n_norm = []
        for t in range(nframes-2):
            if t == 0:
                cellstate_1,lineage_1,cellstate_2,lineage_2 = load2pickles(simulation,t)
            else: 
                cellstate_2,lineage_2 = loadPickle_lite(simulation,t+1)
                
            cellstate_1,cellstate_2 = add_proteins_2(cellstate_1,cellstate_2,lineage_2,t,lambd=lambd)
            cellstate_1 = add_radius_angle_area(cellstate_1) #get r for each cell
            
            n_i,n_prot, r_bins = obtain_convergent_curves(cellstate_1,nbins)
            n_n = (n_i-n_prot)
            n_norm.append(n_n/max(n_i))
            
            cellstate_1 = cellstate_2   
            
        n_norm = np.array(n_norm)
        np.savetxt(path_to_write+'/'+folders[i]+'-r-'+str(repet)+'.gz',n_norm)
    i+=1
'''