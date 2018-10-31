import numpy as np
import Conv_curves_lowram as CLR
import matplotlib.pyplot as plt
import Networking as NWK
class Sector:
    def __init__(self,R_max_t,center,theta,beta):
        self.R = R_max_t
        self.center = center
        if beta < theta:
            beta = beta + 2*np.pi
        self.theta = theta
        self.beta = beta
        self.alpha = beta - theta
        self.x_1 = center[0] + R_max_t*np.cos(theta)
        self.y_1 = center[1] + R_max_t*np.sin(theta)
        self.x_2 = center[0] + R_max_t*np.cos(beta)
        self.y_2 = center[1] + R_max_t*np.sin(beta)
        self.cells = {}
    
    def add_cells_to_sector(self,cellstate):
        for id,cell in cellstate.iteritems():
            if self.check_cell_in_sector(cell) == True:
                self.cells[id] = cell
        self.ncells = len(self.cells)
                
    def check_cell_in_sector(self,cell):
        if self.beta > 2*np.pi:
            beta = self.beta - 2*np.pi
            
            if cell.phi >= self.theta or cell.phi < beta:
                return True
            else:
                return False
        
        if cell.phi >= self.theta  and cell.phi < self.beta:
            return True
        else:
            return False
        
    def resize_sector(self,epsilon,cellstate):
        alpha = epsilon*self.alpha
        delta_alpha = (alpha-self.alpha)
        self.beta = self.beta - delta_alpha/2
        self.theta = self.theta + delta_alpha/2
        self.cells = {}
        self.add_cells_to_sector(cellstate)
        
'''        
fname = "/Users/Medina/cellmodeller/data/Practice_Script_Blank-18-08-21-13-44/step-%05d.pickle"
cellstate,lineage = CLR.loadPickle_lite(fname,1500)
cellstate = CLR.add_radius_angle_area(cellstate)
center = CLR.find_zero(cellstate)
R_max_t = CLR.get_R_max_t(cellstate)

#define beta and theta from a branch
branch_test = NWK.branch_to_cells(cellstate_2,Oak.branch[3])
branch_max = max([cell.r_dist for cell in branch_test])
branch_out = [cell for cell in branch_test if cell.r_dist/branch_max > 0.9]
branch_phi = np.array([cell.phi for cell in branch_out])

avg_phi = np.average(branch_phi)

#caso promedio fuera de la lista

theta = min(branch_phi[branch_phi>avg_phi])
beta = max(branch_phi[branch_phi<avg_phi])

#caso promedio dentro de la lista
#beta = max(branch_phi[branch_phi>avg_phi])
#theta = min(branch_phi[branch_phi<avg_phi])


test_sector = Sector(R_max_t,center,theta,beta)
test_sector.add_cells_to_sector(cellstate)

for cell in branch_test:
    plt.plot(cell.pos[0],cell.pos[1],"ko",markersize = 1)
    
for id,cell in test_sector.cells.iteritems():
    plt.plot(cell.pos[0],cell.pos[1],"bo",markersize = 0.5)
    

epsilon = len(branch_test)/len(test_sector.cells)

'''    
'''
branch_test = NWK.branch_to_cells(cellstate_2,Oak.branch[5])
branch_max = max([cell.r_dist for cell in branch_test])
branch_out = [cell for cell in branch_test if cell.r_dist/branch_max > 0.9]
branch_phi = [cell.phi for cell in branch_out]
beta = min(branch_phi)
theta = max(branch_phi)
'''
