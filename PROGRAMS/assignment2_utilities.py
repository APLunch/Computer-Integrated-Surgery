#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 15:44:57 2021

@author: ychen215
"""

import math
import numpy as np

N = 5

def correction_function(measured_value, coefficient):
    measured_value_scaled = ScaleToBox(measured_value)

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


def bernstein_polynomial(measured_value, ground_true_value, N):    
    
    measured_value_scaled = ScaleToBox(measured_value)

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


def ScaleToBox(q):
    q_min = np.amin(q)
    q_max = np.amax(q)
    
    scaled_q = (q - q_min) / (q_max - q_min)
    
    return scaled_q


def B_Nk(N, k, v):
    N_k = (math.factorial(N))/((math.factorial(k))*(math.factorial(N-k)))
    return N_k * (1-v)**(N-k) * v**k


def F_ijk(N, i, j, k, x, y, z):
    return B_Nk(N, i, x) * B_Nk(N, j, y) * B_Nk(N, k, z)











