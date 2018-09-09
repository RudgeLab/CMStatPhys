import AddProteins as App

from RunAndCompress import GetSubDir


def not_another_main(nframes,lambd):
    
    datafolders = []
    rootdir = "/Users/Medina/cellmodeller/data"
    
    datafolders,datafiles,folders = GetSubDir(rootdir)


    startframe = 0

    
    '''
    foldername = sys.argv[1]
    startframe = sys.argv[2]
    nframes = sys.argv[3]
    dt = sys.argv[4] #There's a bit of trouble with this
    gridsize = sys.argv[5]
    forwards = sys.argv[6]
    '''         
    
    datafolders,datafiles,folders = GetSubDir(rootdir)
    i = 0
    for simulation in datafiles:
        print 'Loading and running '+ datafolders[i]
        cellstates = App.add_protein_pickles(simulation,startframe,nframes,lambd = lambd)
        cellstates_reordered = [cellstates[t]['cellStates'] for t in range(nframes)]
        return cellstates_reordered
    
def load_unmod(nframes):
    
    datafolders = []
    rootdir = "/Users/Medina/cellmodeller/data"
    
    datafolders,datafiles,folders = GetSubDir(rootdir)


    startframe = 0

    
    '''
    foldername = sys.argv[1]
    startframe = sys.argv[2]
    nframes = sys.argv[3]
    dt = sys.argv[4] #There's a bit of trouble with this
    gridsize = sys.argv[5]
    forwards = sys.argv[6]
    '''         
    
    datafolders,datafiles,folders = GetSubDir(rootdir)
    i = 0
    for simulation in datafiles:
        print 'Loading and running '+ datafolders[i]
        cellstate,lineage = App.loadPickle_pro(simulation,startframe,nframes,1,True)
        cellstate = [cellstates[t]['cellStates'] for t in range(nframes)]
        return cellstate