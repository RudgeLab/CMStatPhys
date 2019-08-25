import numpy as np
import matplotlib.pyplot as plt
import os

def load_arrays(path):
    array = np.loadtxt(path)
    return array

t1list = [1300]
#t1 = 1000
#t2 = t1+200
ttree = 1700
root = "/Users/Medina/cellmodeller/histograms_1700"
fig, ax = plt.subplots()
for t1 in t1list:
    t2 = 1700
    fig, ax = plt.subplots()

    Gammapth = root + "/Gamma/Gamma-"+str(t1)+"-"+str(t2)+"-"+str(ttree)
    Gamma = load_arrays(Gammapth)
    
    plt.xlabel("$log_{10}(\Gamma)$")
    plt.ylabel("$log_{10}(freq)$")
    weights = np.ones_like(np.log10(Gamma))/len(np.log10(Gamma))
    ax.hist(np.log10(Gamma), edgecolor = "black",weights=weights,log = True,bins = 10, alpha = 0.5, label='$t_0,t ='+str(t1)+","+str(t2)+'$')

    plt.grid()
    legend = ax.legend(loc='upper right', fontsize=8)

    #plt.hist(np.log10(Gamma), log = True, bins = 10,alpha=0.5)

'''
masterNb = []
for t1 in t1list:
    Nbranchpth = root + "/Nbranch/Nb-"+str(t1)+"-"+str(t2)+"-"+str(ttree)
    Nb = load_arrays(Nbranchpth)
    masterNb.append(Nb)
    plt.hist(np.log10(Nb),bins = 10,alpha=0.5,log = True)
mNb_flat = [item for sublist in masterNb for item in sublist]
'''

'''
masterNn = []

for t1 in t1list:
    Nnormpth = root + "/Nnorm/Nnorm-"+str(t1)+"-"+str(t2)+"-"+str(ttree)
    Nnorm = load_arrays(Nnormpth)
    masterNn.append(Nnorm)
    plt.hist(Nnorm,bins = 10,alpha=0.5,log = True)
mNn_flat = [item for sublist in masterNn for item in sublist]
'''