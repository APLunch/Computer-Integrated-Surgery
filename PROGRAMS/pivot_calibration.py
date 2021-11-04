# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 16:21:15 2021

@author: MikeH
"""
import cismath as cis
import numpy as np
import registration
import assignment2_utilities as pa2
import plotter

#Numpy print options
np.set_printoptions(precision = 2)

def EM_Pivot_Calibration(filename):
    '''
    Pivot calibration for EM Tool

    Parameters
    ----------
    filename : string
        the name of the file that contains data of EM trackers on the probe.

    Returns
    -------
    tuple(Vec3D, Vec3D)
        The position of tool tip and the position of calibration dimple, respectively
    '''
    #Read files and load EM Calibration data data
    NUM_EM_MARKERS = 0
    NUM_EM_DATA_FRAMES = 0
    EM_Data_Frames = []
    
    with open("../Input Data/{}".format(filename),'r') as f:
        data_frame = []
        for i,line in enumerate(f):
            words = line.split()
            #Strip words
            for w in range(len(words)):
                words[w] = words[w].strip(' .,')
            #Handle file header, containing data info
            if i == 0:
                NUM_EM_MARKERS = int(words[0])
                NUM_EM_DATA_FRAMES = int(words[1])
                continue
            else:
                #Handle data
                x, y, z = [float(word) for word in words]
                p = cis.Vec3D(x,y,z)
                data_frame.append(p)
            #Store a frame of data \
            if i % NUM_EM_MARKERS == 0:
                EM_Data_Frames.append(data_frame.copy())
                data_frame.clear()
    
    #Calculate probe position and orientation
    x = sum([point.x for point in EM_Data_Frames[0]])/NUM_EM_MARKERS
    y = sum([point.y for point in EM_Data_Frames[0]])/NUM_EM_MARKERS
    z = sum([point.z for point in EM_Data_Frames[0]])/NUM_EM_MARKERS
    G0 = cis.Vec3D(x,y,z)
    g_list = []
    for col,G_j in enumerate(EM_Data_Frames[0]):
        g_j = G_j-G0
        g_list.append(g_j)
    
    F_list = []
    for data_frame in EM_Data_Frames:
        G_list = []
        for col,G_j in enumerate(data_frame):
            G_list.append(G_j)
        #Find Fk
        F = registration.registration(g_list, G_list)
        F_list.append(F)
        #print(F)
    
    
    #Solve t_G
    A = np.vstack([ np.hstack((F.R.matrix, -1*np.eye(3)))  for F in F_list])
    B = np.vstack([ -1*F.p.matrix for F in F_list])
    X = np.linalg.lstsq(A, B, rcond=None)
    p_pivot = X[0][3:]
    p_tip = X[0][:3]
    
    
    return (cis.Vec3D(p_tip), cis.Vec3D(p_pivot))


def OP_Pivot_Calibration(filename, calbody_filename):
    '''
    Pivot Calibration for Optical Probe

    Parameters
    ----------
    filename : string
        filename of optical tracker data 
    calbody_filename : string
        name of the file containing calbody data.

    Returns
    -------
    tuple(Vec3D, Vec3D)
         Vec3D objects describing the position of tool tip and the calibration post'

    '''
 #Read files and load Optical Calibration data data  
    NUM_OP_MARKERS_BASE = 0
    NUM_OP_MARKERS_PROBE = 0
    NUM_OP_DATA_FRAMES = 0
    OP_Base_Data_Frames = []
    OP_Probe_Data_Frames = []
    OP_Markers_On_Base = []
    
    #Get data header info for 
    with open("../Input Data/{}".format(filename),'r') as f:
        all_lines = f.readlines()
        words = all_lines[0].split()
        words = [word.strip(' ,.') for word in words]
        NUM_OP_MARKERS_BASE = int(words[0])
        NUM_OP_MARKERS_PROBE = int(words[1])
        NUM_OP_DATA_FRAMES = int(words[2])
        
        #Accquire base and probe markers data for each frame
        for n in range(NUM_OP_DATA_FRAMES):
            frame_lines = all_lines[n*(NUM_OP_MARKERS_BASE+NUM_OP_MARKERS_PROBE)+1:(n+1)*(NUM_OP_MARKERS_BASE+NUM_OP_MARKERS_PROBE)+1]
            base_markers_lines = frame_lines[:NUM_OP_MARKERS_BASE]
            #Base markers
            base_frame = []
            for line in base_markers_lines:
                  pos = [ float(word.strip(' ,.')) for word in line.split()]
                  base_marker_pt = cis.Vec3D(pos[0],pos[1],pos[2])
                  base_frame.append(base_marker_pt)
            #Probe markers
            probe_frame = []
            probe_markers_lines = frame_lines[NUM_OP_MARKERS_BASE:]
            for line in probe_markers_lines:
                  pos = [ float(word.strip(' ,.')) for word in line.split()]
                  probe_marker_pt = cis.Vec3D(pos[0],pos[1],pos[2])
                  probe_frame.append(probe_marker_pt)
            #Add data frame to the frame list
            OP_Base_Data_Frames.append(base_frame.copy())
            OP_Probe_Data_Frames.append(probe_frame.copy())
            
    #Get calbody data 
    with open("../Input Data/{}".format(calbody_filename),'r') as f:
        lines = f.readlines()[1:NUM_OP_MARKERS_BASE+1]
        for line in lines:
            words = line.split()
            words = [word.strip(' ,.') for word in words]
            OP_Markers_On_Base.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
    
    #Get registration from optical to EM in first frame
    FD0_inv = registration.registration(OP_Base_Data_Frames[0], OP_Markers_On_Base)
    
    #Calculate probe position and orientation
    x = sum([point.x for point in OP_Base_Data_Frames[0]])/NUM_OP_MARKERS_PROBE
    y = sum([point.y for point in OP_Base_Data_Frames[0]])/NUM_OP_MARKERS_PROBE
    z = sum([point.z for point in OP_Base_Data_Frames[0]])/NUM_OP_MARKERS_PROBE
    G0 = cis.Vec3D(x,y,z)
    #Change respective to EM Base
    G0 = FD0_inv*G0
    g_list = []
    
    for col,G_j in enumerate(OP_Probe_Data_Frames[0]):
        g_j = FD0_inv*G_j-G0
        g_list.append(g_j)
    
    F_list = []
    for n,data_frame in enumerate(OP_Probe_Data_Frames):
        FD_inv = registration.registration(OP_Base_Data_Frames[n], OP_Markers_On_Base)
        G_list = []
        for col,G_j in enumerate(data_frame):
            G_list.append(FD_inv*G_j)
        #Find Fk
        F = registration.registration(g_list, G_list)
        F_list.append(F)
    
    
    #Solve t_G
    A = np.vstack([ np.hstack((F.R.matrix, -1*np.eye(3)))  for F in F_list])
    B = np.vstack([ -1*F.p.matrix for F in F_list])
    X = np.linalg.lstsq(A, B,rcond=None)
    pt = X[0][3:]
    p_pivot = X[0][:3]
    return (cis.Vec3D(pt), cis.Vec3D(p_pivot))

