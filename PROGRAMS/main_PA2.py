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
import registration
import plotter

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
p_tip, p_pivot = pivot_calibration.EM_Pivot_Calibration(filename)

print('pt before correction\n',p_tip, '\n')

#p_tip_corrected_mat= pa2.correction_function(p_tip.matrix, coefficient.T, N, scale_box = scale_box )
p_tip_corrected_mat=p_tip.matrix

print('p_pivot before correction\n',p_pivot, '\n')
print('pt after correction\n',p_tip_corrected_mat, '\n')



#=====================Step 4==================================================
print("Step 4: compute the locations of the fiducials point\n")

filename = 'pa2-debug-a-em-fiducialss.txt'
p_tip_corrected = cis.Vec3D(p_tip_corrected_mat)
fiducials_2_Base = pa2.fiducials_relative_base(filename, p_tip_corrected, coefficient.T, N, scale_box)
for f in fiducials_2_Base:
    print(f)

filename = 'pa2-debug-a-ct-fiducials.txt'

fiducials_2_CT = pa2.fiducials_relative_CT(filename)

print('\n\n')
for f in fiducials_2_CT:
    print(f)
#=====================Step 5==================================================
print("Step 5 Compute the registration frame\n")
F_reg = registration.registration(fiducials_2_Base, fiducials_2_CT)

#DEBUG
plotter.plot_vec_list(fiducials_2_Base)
plotter.plot_vec_list(fiducials_2_CT, color = 'r')

test = [F_reg*p for p in fiducials_2_Base]
plotter.plot_vec_list(fiducials_2_CT, color = 'g')

print("F_reg =")
print(F_reg)
     

#=====================Step 6==================================================
print("Step 6 Compute the pointer tip coordinates with respect to the tracker base\n")       
filename = 'pa2-debug-a-EM-nav.txt'


            
            
            
            
            
            
            
            
            
            
            
            
            
            