import Conv_curves_lowram as CLR
from RunAndCompress import GetSubDir
from pandas import DataFrame
import numpy as np
import os

root = "/Users/Medina/cellmodeller"
datadir = root+"/data"
datafolders,datafiles,folders = GetSubDir(datadir)

i = 0
tf = 1500
world = 250.0
image_size = 1024.0

for simulation in datafiles:
    path_to_write = datadir+"/csv/"+folders[i]
    if not os.path.isdir(path_to_write):
        os.makedirs(path_to_write)
    print 'Loading and running '+ datafolders[i]
    ##here inside iterate t
    for t in range(0,tf):
        
        cellstate, lineage = CLR.loadPickle_lite(simulation,t)
        id_list,x_pixel_list,y_pixel_list,z_pixel_list,length_pixel_list,dir_x_list,dir_y_list,dir_z_list,R_list,G_list,B_list,x_list,y_list,z_list,length_list = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
        #x_pixel_list,y_pixel_list,z_pixel_list,length_pixel_list,dir_x_list,dir_y_list,dir_z_list,bw_list,bh_list,R_list,G_list,B_list,x_list,y_list,z_list,length_list = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    
        for id,cell in cellstate.iteritems():
            position = cell.pos
            color = cell.color
            length = cell.length
            direction = cell.dir
            
            factor = image_size/world
            #calculating and unpacking data for saving
            x_pixel,y_pixel,z_pixel =  image_size - (np.array(position)+(world/2))*factor
            length_pixel = cell.length*factor
            dir_x,dir_y,dir_z = cell.dir[0], cell.dir[1], cell.dir[2]
            R,G,B = cell.color
            x,y,z = position
            length = cell.length
            
            bw =  abs(dir_x)*length_pixel
            bh =  abs(dir_y)*length_pixel
            
            id_list.append(id)
            
            x_pixel_list.append(x_pixel)
            y_pixel_list.append(y_pixel)
            z_pixel_list.append(z_pixel-512.0) #hard fixing
            
            length_pixel_list.append(length_pixel)
            
            dir_x_list.append(dir_x)
            dir_y_list.append(dir_y)
            dir_z_list.append(dir_z)
            
            #bw_list.append(bw)
            #bh_list.append(bh)
            
            R_list.append(R)
            G_list.append(G)
            B_list.append(B)
            
            x_list.append(x)
            y_list.append(y)
            z_list.append(z)
            
            length_list.append(length)
            
        dict_by_id = {"id": id_list,
                      "x_pixel": x_pixel_list,
                          "y_pixel": y_pixel_list,
                          "z_pixel": z_pixel_list,
                          "length_pixel": length_pixel_list,
                          "dir_x": dir_x_list,
                          "dir_y": dir_y_list,
                          "dir_z": dir_z_list,
                          "R": R_list,
                          "G": G_list,
                          "B": B_list,
                          "x_CM": x_list,
                          "y_CM": y_list,
                          "z_CM": z_list,
                          "length_CM": length_list,
                          }
        column = ["id","x_pixel","y_pixel","z_pixel","length_pixel","dir_x","dir_y","dir_z","R","G","B","x_CM","y_CM","z_CM","length_CM"]
        DF = DataFrame.from_dict(dict_by_id)
        DF = DF[column]
        DF.to_csv(path_to_write+"/step-%05d" %(t)+".csv")
        
    i += 1