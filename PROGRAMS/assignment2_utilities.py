#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 15:44:57 2021

@author: ychen215
"""

import math
import numpy as np
import cismath as cis
import registration
import plotter

def correction_function_rowvec(measured_value, coefficient, N, scale_box = None):
    
    if scale_box is None:
        box_min = np.amin(measured_value)
        box_max = np.amax(measured_value)
    else:
        box_min = scale_box[0]
        box_max = scale_box[1]
        
    measured_value_scaled = ScaleToBox(measured_value,box_min,box_max)
    x = measured_value_scaled[:, 0]
    y = measured_value_scaled[:, 1]
    z = measured_value_scaled[:, 2]
    
    p = np.zeros([len(measured_value), 3])
    for row in range(len(measured_value)):
        for i in range(N+1):
            for j in range(N+1):
                for k in range(N+1):
                    #print(coefficient[row, :] * pa2.F_ijk(N, i, j, k, x[row], y[row], z[row]))
                    p[row, :] += coefficient[i * (N+1)**2 + j * (N+1) + k, :] * F_ijk(N, i, j, k, x[row], y[row], z[row])
                    
    return p

def correction_function(measured_value, coefficient, N, scale_box = None):
    '''
    Takes in 3 x N matrix that consists of all measured points(column vectors), the distortion coefficient and
    polynomial distortion degree, return the distortion-corrected matrix 

    Parameters
    ----------
    measured_value : 3 x n numpy array
        measured value with distortion.
    coefficient : 3 x m numpy array
        distortion coefficients
    N : int
        degree of polynomial distortion.

    Returns
    -------
    p :  3 x n numpy array
        matrix of vectors after distorton correction

    '''
    #Convert column vec matrix to row vector matrix
    measured_value_rowvec = measured_value.T
    coefficient_T_rowvec = coefficient.T
    
    p_rowvec = correction_function_rowvec(measured_value_rowvec,coefficient_T_rowvec,N,scale_box = scale_box)

                    
    return p_rowvec.T





def bernstein_polynomial(measured_value, ground_true_value, N, scale_box = None):    
    
    if scale_box is None:
        box_min = np.amin(measured_value)
        box_max = np.amax(measured_value)
    else:
        box_min = scale_box[0]
        box_max = scale_box[1]
        
    measured_value_scaled = ScaleToBox(measured_value,box_min,box_max)

    x = measured_value_scaled[:, 0]
    y = measured_value_scaled[:, 1]
    z = measured_value_scaled[:, 2]
    
    
    A = np.zeros([len(measured_value), (N+1)**3])
    for row in range(len(measured_value)):
        for i in range(N+1):
            for j in range(N+1):
                for k in range(N+1):
                    A[row, i * (N+1)**2 + j * (N+1) + k] = F_ijk(N, i, j, k, x[row], y[row], z[row])
    
    m = np.linalg.lstsq(A, ground_true_value,rcond=None)
    
    coefficient = m[0]
    
    return coefficient


def ScaleToBox(q, q_min, q_max):
    '''
    given a point(vector) q, scale it to the box defined by q_min and q_max

    Parameters
    ----------
    q : 3 x 1 numpy array
        vector to be scaled
    q_min : float
        lower bound
    q_max : float
        upper bound

    Returns
    -------
    scaled_q : 3 x 1 numpy array

    '''
    scaled_q = (q - q_min) / (q_max - q_min)
    return scaled_q




def B_Nk(N, k, v):
    N_k = (math.factorial(N))/((math.factorial(k))*(math.factorial(N-k)))
    return N_k * (1-v)**(N-k) * v**k


def F_ijk(N, i, j, k, x, y, z):
    return B_Nk(N, i, x) * B_Nk(N, j, y) * B_Nk(N, k, z)


def fiducials_relative_base(filename, pt, dist_coefficients, N, scale_box):
    '''
    This function finds the fiducials' locations from a given file and a specified p_tip vector

    Parameters
    ----------
    filename : str
        name of the file that contains fiducial positions relative to CT model
    pt : cismath.Vec3D
        precise displacement of probe tip from probe frame

    Returns
    -------
    fiducial_list : TYPE

    '''
    NUM_EM_MARKERS = 0
    NUM_EM_DATA_FRAMES = 0
    EM_Data_Frames = []
    
    with open("../Input Data/{}".format(filename),'r') as f:
        data_frame = []
        for i,line in enumerate(f):
            words = line.split(',')
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
            #Store a frame of data 
            if i % NUM_EM_MARKERS == 0:
                #Distortion Correction
                '''
                mat = cis.vec_list_to_matrix(data_frame)
                mat_corrected = correction_function(mat, dist_coefficients, N, scale_box = scale_box )
                data_frame = cis.matrix_to_vec_list(mat_corrected)
                '''
                EM_Data_Frames.append(data_frame.copy())
                data_frame.clear()
                
    #Calculate probe frmae and orientation
    x = sum([point.x for point in EM_Data_Frames[0]])/NUM_EM_MARKERS
    y = sum([point.y for point in EM_Data_Frames[0]])/NUM_EM_MARKERS
    z = sum([point.z for point in EM_Data_Frames[0]])/NUM_EM_MARKERS
    G0 = cis.Vec3D(x,y,z)
    g_list = []
    for col,G_j in enumerate(EM_Data_Frames[0]):
        g_j = G_j-G0
        g_list.append(g_j)
        
    #Calculate pose of probe for each probe
    F_list = []
    for data_frame in EM_Data_Frames:
        G_list = []
        for col,G_j in enumerate(data_frame):
            G_list.append(G_j)
        #Find Fk
        F = registration.registration(g_list, G_list)
        F_list.append(F)
        
    
    #Calculate fiducial locations
    fiducial_list = []
    for F in F_list:
        fid_pos = F*pt
        fiducial_list.append(fid_pos)
        
    return fiducial_list


def fiducials_relative_CT(filename):
    with open("../Input Data/{}".format(filename),'r') as f:
            ct_fiducials = []
            for i,line in enumerate(f):
                words = line.split(',')
                #Strip words
                for w in range(len(words)):
                    words[w] = words[w].strip(' .,')
                #Handle file header, containing data info
                if i == 0:
                    NUM_CT_FIDUCIALS = int(words[0])
                    continue
                else:
                    #Handle data
                    x, y, z = [float(word) for word in words]
                    p = cis.Vec3D(x,y,z)
                    ct_fiducials.append(p)
    
    return ct_fiducials



def calc_nav_points(filename, pt, dist_coefficients, N, scale_box, F_reg):
    NUM_EM_MARKERS = 0
    NUM_EM_DATA_FRAMES = 0
    EM_Data_Frames = []
    
    with open("../Input Data/{}".format(filename),'r') as f:
        data_frame = []
        for i,line in enumerate(f):
            words = line.split(',')
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
                #Distortion Correction
                '''
                mat = cis.vec_list_to_matrix(data_frame)
                mat_corrected = correction_function(mat, dist_coefficients, N, scale_box = scale_box )
                data_frame = cis.matrix_to_vec_list(mat_corrected)
                '''
                EM_Data_Frames.append(data_frame.copy())
                data_frame.clear()
                
    #Calculate probe frmae and orientation
    x = sum([point.x for point in EM_Data_Frames[0]])/NUM_EM_MARKERS
    y = sum([point.y for point in EM_Data_Frames[0]])/NUM_EM_MARKERS
    z = sum([point.z for point in EM_Data_Frames[0]])/NUM_EM_MARKERS
    G0 = cis.Vec3D(x,y,z)
    g_list = []
    for col,G_j in enumerate(EM_Data_Frames[0]):
        g_j = G_j-G0
        g_list.append(g_j)
        
    #Calculate pose of probe for each frame
    F_list = []
    for data_frame in EM_Data_Frames:
        G_list = []
        for col,G_j in enumerate(data_frame):
            G_list.append(G_j)
        #Find Fk
        F = registration.registration(g_list, G_list)
        F_list.append(F)
    
    #Calculate p_tip for each frame
    nav_list = []
    for F in F_list:
        nav = F*pt
        nav_list.append(F_reg*nav)
    return nav_list


