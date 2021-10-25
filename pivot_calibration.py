# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 16:21:15 2021

@author: MikeH
"""
import cismath as cis
import numpy as np
import registration
import matplotlib.pyplot as plt
import random

#Numpy print options
np.set_printoptions(precision = 3)

#Read files and load EM Calibration data data
NUM_EM_MARKERS = 0
NUM_EM_DATA_FRAMES = 0
EM_Data_Frames = []

with open("Data/pa1-debug-a-empivot.txt",'r') as f:
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


#Solve t_G
A = np.vstack([ np.hstack((F.R.matrix, -1*np.eye(3)))  for F in F_list])
B = np.vstack([ -1*F.p.matrix for F in F_list])
X = np.linalg.lstsq(A, B)
print(X)