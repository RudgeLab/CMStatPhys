import numpy as np
import matplotlib.pyplot as plt
import os

def load_arrays(path):
    array = np.loadtxt(path)
    return array

t1list = [800,1000,1200,1400]
t1 = 1000
t2 = 1700
ttree = 1700
root = "/Users/Medina/cellmodeller/histograms_1700"

for t1 in t1list:
    Gammapth = root + "/Gamma/Gamma-"+str(t1)+"-"+str(t2)+"-"+str(ttree)
    Gamma = load_arrays(Gammapth)
    plt.hist(np.log10(Gamma),log = True, bins = 15,alpha=0.5)

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