import cPickle
import numpy as np
import matplotlib.pyplot as plt
import Conv_curves_lowram as CLR

def load2pickles(fname,j): #loads pickle j and j+1
    data_1 = cPickle.load(open(fname%(j))) 
    data_2 = cPickle.load(open(fname%(j+1))) 

    cellstate_1 = data_1['cellStates'] 
    lineage_1 = data_1['lineage']
    cellstate_2 = data_2['cellStates'] 
    lineage_2 = data_2['lineage']
    
    return cellstate_1,lineage_1,cellstate_2,lineage_2

def loadPickle_lite(fname,j): #loads pickle j and j+1
    data = cPickle.load(open(fname%(j))) 

    cellstate = data['cellStates'] 
    lineage = data['lineage']

    return cellstate,lineage

def plot_cell_center(cellstate):
    for id,cell in cellstate.iteritems():
        x,y,z = cell.pos
        plt.plot(x,y,'ro',markersize = 0.2)
        
class Create_Tree:
    def __init__(self,lineage):
        self.lineage = lineage
        inv_map = {}
        for did,pid in lineage.iteritems():
            inv_map[pid] = inv_map.get(pid, [])
            inv_map[pid].append(did)
        self.inv_lin = inv_map
        self.branch = {}     
    def create_all_branches(self):
        self.branch[1] = Branch(1)
        self.branch[1].add_nodes_and_edges(self.inv_lin)
        for id,pid in self.lineage.iteritems():
            self.branch[id] = Branch(id)
            self.branch[id].add_nodes_and_edges(self.inv_lin)
            
    def add_t0_branches(self,cellstate,t):
        for bid,branch in self.branch.iteritems():
            if branch.t0 == -1 and cellstate.has_key(bid):
                branch.t0 = t  
    def set_t0_branches(self,cellstate,t):
        for id,cell in cellstate.iteritems():
            if self.branch.has_key(id):
                self.branch[id].t0 = t  
            else:
                print "careful there, branch",id,"doesn't exist in",t
                
    def set_r0_branches(self,cellstate,t):
        for id,cell in cellstate.iteritems():
            if self.branch.has_key(id):
                if self.branch[id].r0 == -1:
                    self.branch[id].r0 = cell.r_dist
            else:
                print "caregul there branch",id,"doesn't exist in",t
            
            
        
class Branch:
    def __init__(self,id):
        self.id = id
        self.nodes = [id]
        self.edges = {}
        self.pos = {}
        self.t0 = -1
        self.r0 = -1
    def add_nodes_and_edges(self,inv_lin):
        for id,dids in inv_lin.iteritems():
            if id in self.nodes:
                self.nodes.append(dids[0])
                self.nodes.append(dids[1])
                self.edges[id] = [id,dids[0],dids[1]]
        
    
        
def convergent_branch(Tree,cellstate,branch_id,nbins = 25.0,color = 'k'):
    branch_cell = {}
    R_max_t = CLR.get_R_max_t(cellstate)
    N_prot = CLR.calculate_sum_prot(cellstate)
    n_prot = N_prot/((np.pi*(R_max_t**2)))
    for id in Tree.branch[branch_id].nodes:
        try:
            branch_cell[id] = cellstate[id]
        except:
            a = 0
    n_i,r_bins = CLR.obtain_convergent_curves_branch(branch_cell,R_max_t,nbins = 25.0)
    
    n_norm = (n_i-n_prot)/max(n_i)
    
    return n_norm

def branch_to_cells(cellstate,branch):
    branch_cells = [cellstate[node] for node in branch.nodes if node in cellstate.keys()]
    return branch_cells     

        
        
        

'''
fname = "/Users/Medina/cellmodeller/data/Practice_Script_Blank-18-08-21-13-44/step-%05d.pickle"
#Practice_Script_Blank-18-08-21-14-28  Practice_Script_Blank-18-08-21-15-12
nframes = 1500
lambd = 1.0
dic_pos = {}
dic_prot = {}
dic_ndiv = {}
for t in range(nframes):
    print 'l-----',t
    if t == 0:
        cellstate_1,lineage_1,cellstate_2,lineage_2 = load2pickles(fname,t)
        
    else:
        lineage = lineage_2
        cellstate_2,lineage_2 = loadPickle_lite(fname,t+1)
    
    #cellstate_1,cellstate_2 = CLR.add_proteins_2(cellstate_1,cellstate_2,lineage_2,t,lambd=lambd)
    cellstate_1 = CLR.add_radius_angle_area(cellstate_1)
    cellstate_2 = CLR.add_radius_angle_area(cellstate_2)
    #cellstate_1,cellstate_2 = CLR.add_ndiv(cellstate_1,cellstate_2)
    
    #create dictionary before dumping data
    for id,cell in cellstate_1.iteritems():
        try: #if it divides the last position is the branch node
            dic_pos[id] = cell.pos
            #dic_prot[id] = cell.red_protein
            #dic_ndiv[id] = cell.ndiv
        except KeyError: #if it doesnt divide the last registered position is the node
            dic_pos[id] = cell.pos
            #dic_prot[id] = cell.red_protein
            #dic_ndiv[id] = cell.ndiv
    #for id,cell in cellstate_2.iteritems(): #if histogram fix is used comment these
        #dic_ndiv[id] = cell.ndiv

    cellstate_1 = cellstate_2

print "creating Tree"
Oak = Create_Tree(lineage_2)
Oak.create_all_branches()
'''

        
'''
plot_cell_center(cellstate_2)
'''

#BRANCH SIZE HISTOGRAM
'''
ak = float(len(cellstate_2))
hist = []
for bid in dic_pos.keys():
    real = 0.0
    expected = 1.0/dic_ndiv[bid]
    try:
        for id in Oak.branch[bid].nodes:
            try:
                cellstate_2[id]
                real += 1.0
            except:
                a = 0
        realvalue = real/ak
        print realvalue, expected
        Gamma = realvalue/expected
        hist.append(Gamma)
    except:
        a = 0
plt.hist(np.log10(hist),bins = 100)
'''
#Branch size histogram fix?
'''
t = 100
cellstate_t,lineage_t = loadPickle_lite(fname,t)

ncells_tf = float(len(cellstate_2))

hist = []
for bid in cellstate_t.iteritems():
    real = 0.0
    expected = 1.0/len(cellstate_t)
    try:
        branch_cells = [cellstate_2[id] for id in Oak.branch[bid].nodes]
        real = float(len(branch_cells))
        realvalue = real/ncells_tf
        print realvalue, expected
        Gamma = realvalue/expected
        hist.append(Gamma)
    except:
        a = 0
plt.hist(np.log10(hist),bins = 100)
'''
#NETWORK PLOTS
'''
for id,pid in lineage_2.iteritems():
    try:
        lineage[id] = [pid,dic_pos[id],dic_pos[pid]]
    except KeyError:
        a = 0
for id,listing in lineage.iteritems():
    pid,pos1,pos2 = listing
    x_1,y_1,z_1 = pos1
    x_2,y_2,z_2 = pos2
    plt.plot([x_1,x_2],[y_1,y_2],"r",linewidth = 0.2)
    
    for ix in [8,9,10,11,12,13,14,15]:
        color = np.random.rand(3,1)
        if id in Oak.branch[ix].edges.keys():
            try:
                x_3,y_3,z_3 = listing[1] 
                x_4,y_4,z_4 = listing[2]
                print 'ak'

                plt.plot([x_3,x_4],[y_3,y_4],color,linewidth = 0.5)
            except:
                a = 0
'''