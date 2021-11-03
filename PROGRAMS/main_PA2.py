#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 16:41:06 2021

@author: ychen215
"""

import assignment1_utilities as pa1
import assignment2_utilities as pa2
import pivot_calibration
#import PA2_2 as pa2 
import numpy as np
import cismath as cis

#=========================Step 1===============================================
print("Step 1: determine the value of C_expected")

calbody_filename = 'pa2-debug-a-calbody.txt'
calreading_filename = 'pa2-debug-a-calreadings.txt'

C_expected, N_C, N_Frames = pa1.compute_C_expected(calbody_filename, calreading_filename)

D_list, A_list, C_list, N_D, N_A, N_C, N_Frames = pa1.read_calreadings(calreading_filename)
C = np.transpose(cis.vec_list_to_matrix(C_list))


#=========================Step 2===============================================
print("Step 2: Produce a suitable distortion correction function")

N = 5

#Get scale box for calculation
scale_box = (np.amin(C_expected), np.amax(C_expected))

#calculate the bernstein polynomial to fit data
coefficient = pa2.bernstein_polynomial(C_expected, C, N, scale_box = scale_box )

# use the correction function
p = pa2.correction_function(C_expected.T, coefficient.T, N, scale_box = scale_box )


#=======================Step 3================================================
print("Step 3: Use the distortion correction function to repeat pivot calibration")
filename = 'pa2-debug-a-empivot.txt'
pt, p_pivot = pivot_calibration.EM_Pivot_Calibration(filename)
print(pt, '\n',p_pivot)

pt2= pa2.correction_function(pt.matrix, coefficient.T, N, scale_box = scale_box )

print(pt2, '\n')