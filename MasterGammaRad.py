import numpy as np
import matplotlib.pyplot as plt
import os

def load_arrays_pro(path,nrbins):
    array = np.loadtxt(path)
    new_array = bin_array(array,nrbins)
    return new_array

def bin_array(array,nrbins):
    array_g = array[0]
    array_r = array[1]
    array_r = array_r/max(array_r)

    
    r_bins = np.linspace(0,1,nrbins)
        
    new_array_g= np.zeros(nrbins)
    
    for j in range(0,len(r_bins)-1):
        
        bin_value = 0
        extra_bin = 0
        n = 0
        e = 0
        
        for k in range(len(array_g)):
                        
            if array_r[k] >= r_bins[j] and array_r[k] < r_bins[j+1]: 
                
                bin_value += array_g[k]
                n += 1
                
            if j == len(r_bins)-2 and array_r[k] >= r_bins[j+1]:

                extra_bin += array_g[k]
                e +=1
                
            else:
                bin_value += 0.0
                
        if n != 0:
            bin_value = bin_value/n
            
        new_array_g[j] = bin_value
        
        if j == len(r_bins)-2:
            new_array_g[j+1] = extra_bin/e
        
    return_array = np.array([new_array_g, r_bins])
    return return_array
        
    
#t1list = [100,200,300,400,500,600,700,800,900,1100,1200,1300,1400,1500,1600]
t1list = [200,300,500,700,900,1100,1300,1500]

t2 = 1700
ttree = 1700
#root = "/Users/Medina/cellmodeller/GammaRad"
root = "/Users/Medina/cellmodeller/ViscosityGamma"

fig, ax = plt.subplots()
nrbins = 20
for t1 in t1list:
    #GRpath= root + "/GammaRad-"+str(t)+"-"+str(t2)+"-"+str(t_tree)+"-"+str(i)

    master_GammaRad = np.array([load_arrays_pro(root + "/GammaRad-"+str(t1)+"-"+str(t2)+"-"+str(ttree)+"-"+str(i),nrbins) for i in range(0,20)])
    
    flattened_master = np.average(master_GammaRad,axis = 0)

    Gamma = flattened_master[0]
    Rad = flattened_master[1]
    
    r = Rad/max(Rad)
    
    plt.ylabel("$\Gamma$")
    plt.xlabel("$R_{max} - R$")
    ax.plot(r,np.log(Gamma), alpha = 0.5, label='$t_0 ='+str(t1)+'$')
    '''
    c = 1- float(t1)/(t2-t1)
    pf = np.polyfit(r,Gamma,12)
    p = np.poly1d(pf)
    pr = p(r)
    ax.plot(r,np.log10(pr))
    '''
    plt.grid()
    legend = ax.legend(loc='upper right', fontsize=8)
    
    
    #plt.hist(np.log10(Gamma), log = True, bins = 10,alpha=0.5)
