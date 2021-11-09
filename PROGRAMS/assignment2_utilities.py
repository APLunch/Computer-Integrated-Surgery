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

def correction_function_vec_list(measured_value, coefficient, N, scale_box = None):
    '''
    Same as correction_function() but takes in list of cismath.Vec3 objects as input and output 

    Parameters
    ----------
    measured_value : list(cismath.Vec3)
        measured value with distortion.
    coefficient : 3 x m numpy array
        distortion coefficients
    N : TYPE
        degree of polynomial distortion.
    scale_box : (float, float), optional
        Scale box of bertein polynomial function. The default is None.

    Returns
    -------
    list(cismath.Vec3)
     corrected vectors

    '''
    measured_value_mat = cis.vec_list_to_matrix(measured_value)
    p_mat = correction_function(measured_value_mat, coefficient, N, scale_box = scale_box)
    return cis.matrix_to_vec_list(p_mat)


# calculate the coefficient of bernstein polynomial used in the correctio function
def bernstein_polynomial(ground_true_value, measured_value, N, scale_box = None):    
    
    #scale the input value to the range between 0 and 1
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
    
    # use least square error to find the coefficient
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



#Bernstein basis polynomials of degree N
def B_Nk(N, k, v):
    N_k = (math.factorial(N))/((math.factorial(k))*(math.factorial(N-k)))
    return N_k * (1-v)**(N-k) * v**k

# "tensor" form polynomial using Nth degree Bernstein polynomial
def F_ijk(N, i, j, k, x, y, z):
    return B_Nk(N, i, x) * B_Nk(N, j, y) * B_Nk(N, k, z)


def fiducials_relative_base(filename, calib_local_frame, pt, dist_coefficients, N, scale_box, correction = True):
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
                if correction:
                    #Distortion Correction
                    data_frame = correction_function_vec_list(data_frame, dist_coefficients, N, scale_box = scale_box)
                EM_Data_Frames.append(data_frame.copy())
                data_frame.clear()
                
    #Get local frame
    g_list = calib_local_frame
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



def calc_nav_points(filename,calib_local_frame ,pt, dist_coefficients, N, scale_box, F_reg, correction = True):
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
                if correction:
                    #Distortion Correction
                    data_frame = correction_function_vec_list(data_frame, dist_coefficients, N, scale_box = scale_box)
                EM_Data_Frames.append(data_frame.copy())
                data_frame.clear()
                
    #Get local frame
    g_list = calib_local_frame
        
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
        nav_list.append(F_reg*F*pt)
    return nav_list

# read the calreadings txt file
def read_calreadings(filename):
    
    with open("../Input Data/{}".format(filename),'r') as f:
        lines = f.readlines()
        
        #read first line to get ND, NA and NC
        first_line = lines[0].split(",")
        
        N_D = int(first_line[0])
        N_A = int(first_line[1])
        N_C = int(first_line[2])
        N_Frames = int(first_line[3])
        
        #save all the coordinates of D, A, C to three defined lists of Vec3D object
        D_list = []
        A_list = []
        C_list = []
        
        for i in range(N_Frames):
            for line in lines[1+i*(N_D+N_A+N_C) : N_D+1+i*(N_D+N_A+N_C)]:
                words = line.split()
                words = [word.strip(',') for word in words]
                D_list.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
        
            for line in lines[1+N_D+i*(N_D+N_A+N_C) : N_D+1+N_A+i*(N_D+N_A+N_C)]:
                words = line.split()
                words = [word.strip(',') for word in words]
                A_list.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
        
            for line in lines[N_D+1+N_A+i*(N_D+N_A+N_C) : N_D+1+N_A+N_C+i*(N_D+N_A+N_C)]:
                words = line.split()
                words = [word.strip(',') for word in words]
                C_list.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
        
    return D_list, A_list, C_list, N_D, N_A, N_C, N_Frames



# This function compute the C_expected value as required in Question4
def compute_C_expected(filename1, filename2):  
    # read the coordinates of a and A 
    d_list, a_list, c_list = read_calbody(filename1)
    D_list, A_list, C_list, N_D, N_A, N_C, N_Frames = read_calreadings(filename2)
    
    # initialize the error in Fa and Fd for unit test
    error_Fa_x = 0
    error_Fa_y = 0
    error_Fa_z = 0
    error_Fd_x = 0
    error_Fd_y = 0
    error_Fd_z = 0

    
    #calculate the C_expected
    C_expected_list = []
    for i in range(N_Frames):
        # calculate the corresponding Fd, Fa, and Fc in each frame
        D_sublist = D_list[i*N_D:(i+1)*N_D]
        A_sublist = A_list[i*N_A:(i+1)*N_A]
        C_sublist = C_list[i*N_C:(i+1)*N_C]
        Fd = registration.registration(d_list, D_sublist)
        Fa = registration.registration(a_list, A_sublist)
        Fc = registration.registration(c_list, C_sublist)
        
    
        # calculate the C_expected by using the corresponding Fd, Fa
        for c_vector in c_list:
            c = Fd.inv()*(Fa * c_vector)
            C_expected_list.append(c)
    
        # Unit test for Fa and Fd by calculating the error in Fa and Fd
        A_test = []
        D_test = []
    
        # calcule late error in Fa
        for n in range(len(a_list)):
            a_vector = a_list[n]
            A_test.append(Fa * a_vector)
            error_Fa_x += abs(A_sublist[n].x - (Fa * a_vector).x)
            error_Fa_y += abs(A_sublist[n].y - (Fa * a_vector).y)
            error_Fa_z += abs(A_sublist[n].z - (Fa * a_vector).z)
            
            
        # calcule late error in Fd
        for n in range(len(d_list)):
            d_vector = d_list[n]
            D_test.append(Fd * d_vector)
            error_Fd_x += abs(D_sublist[n].x - (Fd * d_vector).x)
            error_Fd_y += abs(D_sublist[n].y - (Fd * d_vector).y)
            error_Fd_z += abs(D_sublist[n].z - (Fd * d_vector).z)
     
            
    error_Fa_x = error_Fa_x/len(A_list)
    error_Fa_y = error_Fa_y/len(A_list)
    error_Fa_z = error_Fa_z/len(A_list)
    #print('error Fa = ', error_Fa_x, error_Fa_y, error_Fa_z)
    
    
    error_Fd_x = error_Fd_x/len(D_list)
    error_Fd_y = error_Fd_y/len(D_list)
    error_Fd_z = error_Fd_z/len(D_list)
    #print('error Fd = ', error_Fd_x, error_Fd_y, error_Fd_z)
    print('\n')    
    
    
    print("C_expected:")    
    # Save C_expected as an numpy.ndarray
    C_expected = np.transpose(cis.vec_list_to_matrix(C_expected_list))
    
    print(C_expected)
    print('\n')
    
    return C_expected, N_C, N_Frames

# read the calbody txt file
def read_calbody(filename):
    
    with open("../Input Data/{}".format(filename),'r') as f:
        lines = f.readlines()
        
        #read first line to get ND, NA and NC
        first_line = lines[0].split(",")
        N_D = int(first_line[0])
        N_A = int(first_line[1])
        N_C = int(first_line[2])
        
        #save all the coordinates of d, a, c to three defined lists of Vec3D object
        points_d = []
        points_a = []
        points_c = []
        
        for line in lines[1:N_D+1]:
            words = line.split()
            words = [word.strip(',') for word in words]
            points_d.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
    
        for line in lines[1+N_D:N_D+1+N_A]:
            words = line.split()
            words = [word.strip(',') for word in words]
            points_a.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
    
        for line in lines[N_D+1+N_A:N_D+1+N_A+N_C]:
            words = line.split()
            words = [word.strip(',') for word in words]
            points_c.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
    
    return points_d, points_a, points_c