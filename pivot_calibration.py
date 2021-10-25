# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import cismath as cis
import numpy as np

#Numpy print options
np.set_printoptions(precision = 3)

#Read files and load EM Calibration data data
NUM_EM_MARKERS = 0
NUM_EM_DATA_FRAMES = 0
EM_Data_Frames = []

with open("Data/pa1-unknown-h-empivot.txt",'r') as f:
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
F_list = []
for data_frame in EM_Data_Frames:
    x = sum([point.x for point in data_frame])/NUM_EM_MARKERS
    y = sum([point.y for point in data_frame])/NUM_EM_MARKERS
    z = sum([point.z for point in data_frame])/NUM_EM_MARKERS
    G0 = cis.Vec3D(x,y,z)
    g_matrix = np.ones([4,NUM_EM_MARKERS])
    G_matrix = np.ones([4,NUM_EM_MARKERS])
    for col,G_j in enumerate(data_frame):
        g_j = G_j-G0
        g_matrix[:3,col] =  g_j.matrix.reshape(3)
        G_matrix[:3,col] =  G_j.matrix.reshape(3)
    #Find Fk
    F = np.linalg.lstsq(g_matrix, G_matrix)
    print(F)
    