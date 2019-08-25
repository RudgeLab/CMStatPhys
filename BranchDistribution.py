import Networking as NWK
import numpy as np
import matplotlib.pyplot as plt
#from scipy.signal import savgol_filter


def mean_sq_d(array2d):
    if len(array2d.shape) == 2:
        x = array2d[:,0]
        y = array2d[:,1]
        n = array2d.shape[0]
        avg_x = np.average(x)
        avg_y = np.average(y)
        x_dev = x-avg_x
        y_dev = y-avg_y
        total_sum = np.sum(x_dev**2 + y_dev**2)
        MSD = total_sum/n
        return MSD
    else:
        return 0.0

fname = "/Users/Medina/cellmodeller/data/Practice_Script_Blank-18-08-21-13-44/step-%05d.pickle"
t1 = 255
t2 = 1000

cellstate_0,lin_0 = NWK.loadPickle_lite(fname,t1)
bnumber = len(cellstate_0)
Oak.set_t0_branches(cellstate_0,t1)
bid_array = [bid for bid,branch in Oak.branch.iteritems() if branch.t0 == t1]

variances = np.zeros((bnumber,t2-t1))
ncells = np.zeros((bnumber,t2-t1))

#need to find new t0 with n cells which i will track the branches

for t in range(t1,t2):
    print 'v----',t
    cellstate,lineage = NWK.loadPickle_lite(fname,t)    
    i = 0
    for bid in bid_array:
    
        branch = Oak.branch[bid]
        x_t_list = []
        ncell = 0
        for node in branch.nodes:
            if cellstate.has_key(node):
                x = cellstate[node].pos[0]
                y = cellstate[node].pos[1]
                x_t_list.append([x,y])
                ncell += 1
        x_t_array = np.array(x_t_list)
        var_t = mean_sq_d(x_t_array)
        variances[i,t-t1] = var_t
        ncells[i,t-t1] = ncell
        i +=1



print "plotting"
'''
for bid in range(0,bnumber-1):
    color = (np.random.rand(1)[0],np.random.rand(1)[0],np.random.rand(1)[0])
    
    variances_branch = variances[bid+1]
    
    times = np.arange(t1,t2)
    log_time = np.log10(times)
    log_var = np.log10(variances_branch)
    plt.plot(log_time,log_var,color = color,linewidth = 1)
    
    times_2 = times[variances_branch>0]
    variances_branch_2 = variances_branch[variances_branch>0]
    if len(variances_branch_2)>20:
        polyfit = np.polyfit(np.log10(times_2),np.log10(variances_branch_2),1)
        slope = polyfit[0]
        c = polyfit[1]
        plt.plot(np.log10(times_2),slope*np.log10(times_2)+c,color=color)
        
    plt.plot(log_time,log_var,linewidth = 1)

plt.plot(log_time,log_time,"k")
'''


for bid in range(0,bnumber-1):
    color = (np.random.rand(1)[0],np.random.rand(1)[0],np.random.rand(1)[0])
    
    variances_branch = variances[bid+1]
    ncells_b = ncells[bid+1]
    times = np.arange(t1,t2)
    log_time = np.log10(times)
    log_var = np.log10(variances_branch)
    
    
    times_2 = times[ncells_b >10]
    variances_branch_2 = variances_branch[ncells_b>10]
    if len(variances_branch_2)>10:
        polyfit = np.polyfit(np.log10(times_2),np.log10(variances_branch_2),1)
        slope = polyfit[0]
        c = polyfit[1]
        plt.plot(np.log10(times_2),slope*np.log10(times_2)+c,color=color)
        
    plt.plot(log_time,log_var,color = color,linewidth = 1)

#plt.plot(log_time,log_time,"k")

#SLOPE OF LOG SPACE  = ALPHA 
#PLOT AVERAGE ALPHA VS POSITION OF CELL_I

slopes = np.zeros(bnumber)
pos_i = np.zeros(bnumber)
cellstate,lineage =  NWK.loadPickle_lite(fname,t1)
cellstate = CLR.add_radius_angle_area(cellstate)
R_max_t = CLR.get_R_max_t(cellstate)
i=0
for bid in bid_array:
    
    cell_i = cellstate[bid]
    variances_branch = variances[i]
    times = np.arange(t1,t2)
    times_2 = times[ncells_b >10]
    variances_branch_2 = variances_branch[ncells_b >10]
    if len(variances_branch_2)>10:
        polyfit = np.polyfit(np.log10(times_2),np.log10(variances_branch_2),1)
        slope = polyfit[0]
        #DY = np.diff(variances_branch)
        #DX = np.diff(times)
        #slope = np.average(DY/DX)
        pos_i[i] = cell_i.r_dist/R_max_t
        slopes[i] = slope
    i+=1



def bin_slopes(pos_i,slopes,nbins):
    r_bins = np.linspace(0,1,nbins)
        
    alpha_i= np.zeros(nbins)
    
    for j in range(0,len(r_bins)-1):
        
        bin_value = 0
        extra_bin = 0
        n = 0
        e = 0
        
        for k in range(len(slopes)):
                        
            if pos_i[k] >= r_bins[j] and pos_i[k] < r_bins[j+1]: 
                
                bin_value += slopes[k]
                n += 1
                
            if j == len(r_bins)-2 and pos_i[k] >= r_bins[j+1]:
                
                extra_bin += slopes[k]
                e +=1
                
            else:
                bin_value += 0.0
                
        if n != 0:
            bin_value = bin_value/n
            
        alpha_i[j] = bin_value
        if j == len(r_bins)-2:
            alpha_i[j+1] = extra_bin/e
        
    return r_bins,alpha_i
        
   
#variance plot
'''
plt.plot(np.gradient(variances[0]))
yhat = savgol_filter(np.gradient(variances[0]), 51, 3)
plt.plot(yhat)
'''


    



