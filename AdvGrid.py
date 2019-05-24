import numpy as np
import cPickle
import infotheory as IT

#import matplotlib.pyplot as plt

class Grid():
    def __init__(self,gx,gy,dgx,dgy,t,dt,center):
        self.grid = []
        self.gx = gx
        self.gy = gy
        self.dgx = dgx
        self.dgy = dgy
        self.t = t
        self.dt = dt
        self.c = center
        for ix in range(gx):
            for iy in range(gy):
                self.grid.append(Ensemble(self.t,ix,iy,ix*dgx-center,iy*dgy-center))
            
    def __getitem__(self,i):

        if len(i) == 2 and type(i[0]) == int:
            ix,iy = i
            return self.grid[ix*self.gx+iy]
        

    def update(self,cellstate_t1,lineage_t1,cellstate_t2,lineage_t2,grid_t2):
        #print '-'*16,'Step ',t,'-'*16
        #Add cells to gridcells :
        self.add_cells_to_ensembles(cellstate_t1)
        #calc velocity of grid at t = 0 to t+1: 
        Total,counted_total,grid_t2 = self.calc_velocity_of_ensembles(cellstate_t1,lineage_t1,cellstate_t2,lineage_t2,grid_t2)
        if Total > counted_total:
            print 'Cells in the void between ensembles: ', Total - counted_total
        return grid_t2
        
    def add_cells_to_ensembles(self,cellstate):
        for ix in range(self.gx):
            for iy in range(self.gy):
                for (id,cell) in cellstate.iteritems():
                    if self[ix,iy].CheckCellInEnsemble(cell,self.dgx,self.dgy) == True:
                        self[ix,iy].addCell(cell,id)
                        
    def calc_velocity_of_ensembles(self,cellstate_t1,lineage_t1,cellstate_t2,lineage_t2,grid_t2):
        counted_total = 0
        for ix in range(self.gx):
            for iy in range(self.gy):

                Total,cell_no = self[ix,iy].CalcVel(cellstate_t1,cellstate_t2,lineage_t1,self.dt,self.dgx,self.dgy)
                
                counted_total += cell_no
                #Should we remove skipped cells? might affect entropy flux calculation  
                #Move grids t+1:
                grid_t2[ix,iy].px = self[ix,iy].px + self[ix,iy].vx*self.dt
                grid_t2[ix,iy].py = self[ix,iy].py + self[ix,iy].vy*self.dt
                grid_t2[ix,iy].vx = -self[ix,iy].vx
                grid_t2[ix,iy].vy = -self[ix,iy].vy
        return Total, counted_total,grid_t2
    



class Ensemble():
    def __init__(self,t,ix,iy,px0,py0):
        self.px = px0
        self.py = py0
        self.vx = 0
        self.vy = 0
        self.ix = ix
        self.iy = iy
        self.t = t
        self.cells = {}
        self.skipped = 0
        self.cell_number = 0
        self.entropy = {}
        self.averages = {}

    def addCell(self,cell,id): #cell = cellstate
        #print 'added cell: ',id, 'to', (self.px,self.py)
        self.cells[id] = cell
        
    def CalcVel(self,cellstate_t1,cellstate_t2,lineage_t1,dt,dgx,dgy):
        dx,dy = 0,0
        total_cells = 0
        for id,next_cell in cellstate_t2.iteritems():
            dx_cell = 0
            dy_cell = 0
            total_cells += 1
            if self.CheckCellInEnsemble(next_cell,dgx,dgy) == True:
                
                #print 'Calculating velocity of Ensemble: ', (self.px, self.py),' cell id ',id
                try:
                    dx_cell = next_cell.pos[0]-cellstate_t1[id].pos[0]
                    dy_cell = next_cell.pos[1]-cellstate_t1[id].pos[1]
                    self.cell_number += 1
                    #print '+++++Success'
                except KeyError:
                    #print '----- Cell not in current step, checking for division'
                    # Previous cell does not exist, use parent cell
                   # print lineage
                    pids = [key for key, value in lineage_t1.iteritems() if value == id]
                    pid = pids[0]
                    #print 'Using daughter cell: ',pid
                    dx_cell = next_cell.pos[0]-cellstate_t1[pid].pos[0]
                    dy_cell = next_cell.pos[1]-cellstate_t1[pid].pos[1]
                    self.cell_number += 1 # Count as 1/2 to take average of children
            dx += dx_cell
            dy += dy_cell
                        
        if self.cell_number != 0:
            dx = dx/self.cell_number
            dy = dy/self.cell_number
        self.averages['vel'] = [dx/dt,dy/dt]
        self.vx = dx/dt
        self.vy = dy/dt
        return total_cells,self.cell_number
    
    def CheckCellInEnsemble(self,cellstate,dgx,dgy):
        x,y = cellstate.pos[0],cellstate.pos[1]
        xg,yg = self.px, self.py
        if x >= xg and x < xg+dgx and y >= yg and y < yg+dgy:
            return True
        else:
            return False
        
    def calculate_average(self,attribute):
        avg = 0
        n = 0
        for id in self.cells.keys():
            cell_atr = getattr(self.cells[id],attribute, None)
            if cell_atr:
                avg += cell_atr
                n += 1
        self.averages[attribute] = avg/len(self.cells)
'''
def fname2pickle(fname):
    if fname.endswith(".png") or fname.endswith(".jpg"):
        newfname = fname[:len(fname)-4]+".pickle"
    elif fname.endswith(".tiff"):
        newfname = fname[:len(fname)-5]+".pickle"
    return newfname
'''

def loadPickle_lite(fname,j): #loads pickle j and j+1
    data = cPickle.load(open(fname%(j))) 

    cellstate = data['cellStates'] 
    lineage = data['lineage']

    return cellstate,lineage

def timestep(fname,t1,t2,dt,ngrid,size = 100):
    
    cellstate_t1, lineage_t1 = loadPickle_lite(fname,t1)
    cellstate_t2, lineage_t2 = loadPickle_lite(fname,t2)



    gx,gy = ngrid,ngrid
    dgx,dgy = float(size)/ngrid,float(size)/ngrid
    center = size/2
    #print "worldsize = ", worldsize
    #print "resizing factor = ", resizing  
    print "Grid dimensions: ",gx,gy
    
    grid_t1 = Grid(gx,gy,dgx,dgy,t1,dt,center)
    grid_t2 = Grid(gx,gy,dgx,dgy,t2,dt,center)
    
    grid_t2 = grid_t1.update(cellstate_t1,lineage_t1,cellstate_t2,lineage_t2,grid_t2)
    

    return grid_t1,grid_t2,cellstate_t1,cellstate_t2,lineage_t1,lineage_t2

fname = "/Users/Medina/cellmodeller/data/Practice_Script_Blank-18-08-21-13-44/step-%05d.pickle"
t1 = 301
t2 = 300
dt = 1 #There's a bit of trouble with this
ngrid = 10 #pixels per grid
size = 100
agrid_t1,agrid_t2,acellstate_t1,acellstate_t2,alineage_t1,alineage_t2 = timestep(fname,t1,t2,dt,ngrid,size = size)
