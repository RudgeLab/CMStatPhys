import Networking as NWK
import numpy as np
import matplotlib.pyplot as plt
import Conv_curves_lowram as CLR
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
def bin_radius(variances,R_max_t,bid_array,Tree,nbins):
    
    r_bins = np.linspace(0,1,nbins)
    lens = len(variances.values()[0])
    alpha_i= np.zeros((nbins,lens))
    
    for j in range(0,len(r_bins)-1):
        
        bin_value = np.zeros(len(variances.values()[0]))
        extra_bin = np.zeros(len(variances.values()[0]))
        n = 0
        e = 0
        
        for bid,variance_branch in variances.iteritems():
            branch_r0 = Tree.branch[bid].r0/R_max_t
            if branch_r0 >= r_bins[j] and branch_r0 < r_bins[j+1]: 
                
                bin_value += variance_branch
                n += 1
                
            if j == len(r_bins)-2 and branch_r0 >= r_bins[j+1]:
                
                extra_bin += variance_branch
                e +=1
                
            else:
                bin_value += 0.0
                
        if n != 0:
            bin_value = bin_value/n
            
        alpha_i[j] = bin_value
        if j == len(r_bins)-2:
            alpha_i[j+1] = extra_bin/e
        
    return r_bins,alpha_i

def run_script(fname,t1,t2,tm,Tree,nbins):
    
    
    cellstate_0,lin_0 = NWK.loadPickle_lite(fname,t1)
    cellstate_0 = CLR.add_radius_angle_area(cellstate_0)
    Tree.set_t0_branches(cellstate_0,t1)
    Tree.set_r0_branches(cellstate_0,t1)
    bid_array = np.array([bid for bid,branch in Tree.branch.iteritems() if branch.t0 == t1])
    R_max_t = CLR.get_R_max_t(cellstate_0)
    
    variances = {}
    ncells = {}
    for bid in bid_array:
        variances[bid] = np.array([])
        ncells[bid] = np.array([])
    
    #need to find new t0 with n cells which i will track the branches
    
    for t in range(t1,t2):
        print 'v----',t
        cellstate,lineage = NWK.loadPickle_lite(fname,t)    
    
        for bid in bid_array:
            branch = Tree.branch[bid]
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
            variances[bid] = np.append(variances[bid],var_t)
            ncells[bid] = np.append(ncells[bid],ncell)

    r,MSDS = bin_radius(variances,R_max_t,bid_array,Tree,nbins)



    print "plotting"
    
    #this
    alphas = []
    for variance_r in MSDS:
        
        #color = (np.random.rand(1)[0],np.random.rand(1)[0],np.random.rand(1)[0])
        times = np.arange(t1,t2)
        #log_time = np.log10(times)
        #log_var = np.log10(variance_r)
        
        times_2 = times[times>tm]
        variances_r_2 = variance_r[times>tm]
        
        if len(variances_r_2)>10:
            polyfit = np.polyfit(np.log10(times_2),np.log10(variances_r_2),1)
            slope = polyfit[0]
            #c = polyfit[1]
            alphas.append(slope)
            #plt.plot(np.log10(times_2),slope*np.log10(times_2)+c,color=color,linewidth = 0.5)
    return r,MSDS,times,alphas
        #plt.plot(log_time,log_var,color = color,linewidth = 1)
 
#plt.plot(log_time,log_time,"k")
#SLOPE OF LOG SPACE  = ALPHA 
#PLOT AVERAGE ALPHA VS POSITION OF CELL_I
'''
fname = "/Users/Medina/cellmodeller/data/Practice_Script_Blank-18-08-21-13-44/step-%05d.pickle"
t1 = 255
tm = 400
t2 = 1500
nbins = 15
r,MSDS,times,alphas = run_script(fname,t1,t2,tm,Oak,nbins)


slopes = np.zeros(bnumber)
pos_i = np.zeros(bnumber)
cellstate,lineage =  NWK.loadPickle_lite(fname,t1)
cellstate = CLR.add_radius_angle_area(cellstate)
R_max_t = CLR.get_R_max_t(cellstate)
#i=0
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
    #i+=1
'''


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


    



