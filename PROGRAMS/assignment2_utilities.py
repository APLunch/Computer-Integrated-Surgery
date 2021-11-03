#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 15:44:57 2021

@author: ychen215
"""

import math
import numpy as np
import cismath as cis


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
    
    m = np.linalg.lstsq(A, ground_true_value)
    
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











