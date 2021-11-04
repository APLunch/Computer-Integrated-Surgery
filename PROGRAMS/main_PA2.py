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


dataset = 'a'

if dataset in 'abcdefg':
        calbody_filename = 'pa2-debug-'+ dataset +'-calbody.txt'
        calreading_filename = 'pa2-debug-'+ dataset +'-calreadings.txt'
        empivot_filename = 'pa2-debug-' + dataset + '-empivot.txt'
        em_fiducials_filename = 'pa2-debug-' + dataset +'-em-fiducialss.txt'
        ct_fiducials_filename = 'pa2-debug-' + dataset +'-ct-fiducials.txt'
        em_nav_filename = 'pa2-debug-' + dataset +'-EM-nav.txt'
    
elif dataset in 'hijk':
    calbody_filename = 'pa2-unknown-'+ dataset +'-calbody.txt'
    calreading_filename = 'pa2-unknown-'+ dataset +'-calreadings.txt'
    empivot_filename = 'pa2-unknown-' + dataset + '-empivot.txt'
    em_fiducials_filename = 'pa2-unknown-' + dataset +'-em-fiducialss.txt'
    ct_fiducials_filename = 'pa2-unknown-' + dataset +'-ct-fiducials.txt'
    em_nav_filename = 'pa2-unknown-' + dataset +'-EM-nav.txt'



#=========================Step 1===============================================
print("Step 1: determine the value of C_expected")


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
p_tip, p_pivot = pivot_calibration.EM_Pivot_Calibration(empivot_filename)

print('pt before correction\n',p_tip, '\n')

#p_tip_corrected_mat= pa2.correction_function(p_tip.matrix, coefficient.T, N, scale_box = scale_box )
#DEBUG
p_tip_corrected_mat=p_tip.matrix

print('p_pivot before correction\n',p_pivot, '\n')
print('pt after correction\n',p_tip_corrected_mat, '\n')



#=====================Step 4==================================================
print("Step 4: compute the locations of the fiducials point\n")
p_tip_corrected = cis.Vec3D(p_tip_corrected_mat)
fiducials_2_Base = pa2.fiducials_relative_base(em_fiducials_filename, p_tip_corrected, coefficient.T, N, scale_box)
for f in fiducials_2_Base:
    print(f)



fiducials_2_CT = pa2.fiducials_relative_CT(ct_fiducials_filename)

print('\n\n')
for f in fiducials_2_CT:
    print(f)
#=====================Step 5==================================================
print("Step 5 Compute the registration frame\n")
F_reg = registration.registration(fiducials_2_Base, fiducials_2_CT)


#DEBUG
plotter.plot_vec_list(fiducials_2_Base)
plotter.plot_vec_list(fiducials_2_CT, color = 'r')


for i in range(len(fiducials_2_Base)):
    plotter.plot_arrow(fiducials_2_Base[i],fiducials_2_CT[i]-fiducials_2_Base[i])
    

print("F_reg =")
print(F_reg)
     

#=====================Step 6==================================================
print("Step 6 Compute the pointer tip coordinates with respect to the tracker base\n")       
nav_list = pa2.calc_nav_points(em_nav_filename, p_tip_corrected, coefficient.T, N, scale_box, F_reg)

for nav in nav_list:
    print(nav)

            
            
            
            
            
            
            
            
            
            
            
            
            
            