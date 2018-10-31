import Sectors as SCS
import numpy as np
import Conv_curves_lowram as CLR
import matplotlib.pyplot as plt
        
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
    
testangle = [[0,np.pi/2],[0,np.pi/3],[0,np.pi/4],[0,np.pi/8],[0,np.pi],[0,2*np.pi],[0,3*np.pi/4],[np.pi/3,np.pi/2],[np.pi/4,np.pi/3],[np.pi/5,np.pi/4],[np.pi/9,np.pi/8],[np.pi/2,np.pi],[np.pi,2*np.pi],[np.pi/2,3*np.pi/4]]    
for i in range(len(testangle)):
    theta = testangle[i][0]
    beta =  testangle[i][1]
    MeanSqD = np.zeros(1500)
    for t in range(1, 800):
        print '----',t
        fname = "/Users/Medina/cellmodeller/data/Practice_Script_Blank-18-08-21-13-44/step-%05d.pickle"
        cellstate,lineage = CLR.loadPickle_lite(fname,t)
        cellstate = CLR.add_radius_angle_area(cellstate)
        center = CLR.find_zero(cellstate)
        R_max_t = CLR.get_R_max_t(cellstate)
        AlphaSector = SCS.Sector(R_max_t,center,theta,beta)
        AlphaSector.add_cells_to_sector(cellstate)
        cellpos = np.array([[cell.pos[0],cell.pos[1]] for id,cell in AlphaSector.cells.iteritems()])
        MeanSqD[t] = mean_sq_d(cellpos)
    plt.plot(np.log10(np.arange(0,1500)),np.log10(MeanSqD))
    