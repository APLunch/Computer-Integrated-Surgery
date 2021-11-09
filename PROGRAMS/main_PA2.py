#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 16:41:06 2021

@author: ychen215
"""

import assignment2_utilities as pa2
import pivot_calibration
import numpy as np
import cismath as cis
import registration
#import plotter


def main(dataset):
    print('===============================================================')
    print('=======================Data Set={}============================='.format(dataset))
    print('===============================================================')
    
    #Data pre-processing and configuration for filenames
    if dataset in 'abcdef':
            calbody_filename = 'pa2-debug-'+ dataset +'-calbody.txt'
            calreading_filename = 'pa2-debug-'+ dataset +'-calreadings.txt'
            empivot_filename = 'pa2-debug-' + dataset + '-empivot.txt'
            optpivot_filename = 'pa2-debug-' + dataset + '-optpivot.txt'
            em_fiducials_filename = 'pa2-debug-' + dataset +'-em-fiducialss.txt'
            ct_fiducials_filename = 'pa2-debug-' + dataset +'-ct-fiducials.txt'
            em_nav_filename = 'pa2-debug-' + dataset +'-EM-nav.txt'
            output1_filename = 'pa2-debug-'+ dataset +'-output1.txt'
            output2_filename = 'pa2-debug-'+ dataset +'-output2.txt'
            
    elif dataset in 'ghijk':
        calbody_filename = 'pa2-unknown-'+ dataset +'-calbody.txt'
        calreading_filename = 'pa2-unknown-'+ dataset +'-calreadings.txt'
        empivot_filename = 'pa2-unknown-' + dataset + '-empivot.txt'
        optpivot_filename = 'pa2-unknown-' + dataset + '-optpivot.txt'
        em_fiducials_filename = 'pa2-unknown-' + dataset +'-em-fiducialss.txt'
        ct_fiducials_filename = 'pa2-unknown-' + dataset +'-ct-fiducials.txt'
        em_nav_filename = 'pa2-unknown-' + dataset +'-EM-nav.txt'
        output1_filename = 'pa2-unknown-'+ dataset +'-output1.txt'
        output2_filename = 'pa2-unknown-'+ dataset +'-output2.txt'
    
    
    #Loocking for max and minimum value in entire EM dataset, for Berstein polynomial
    p_min = 10000
    p_max = -10000
    for fname in [calreading_filename,empivot_filename,em_fiducials_filename,em_nav_filename]:
        with open("../Input Data/{}".format(fname),'r') as f:
            for i,line in enumerate(f):
                words = line.split(',')
                #Strip words
                for w in range(len(words)):
                    words[w] = words[w].strip(' .,')
                #Handle file header, containing data info
                if i == 0:
                    continue
                else:
                    #Handle data
                    x, y, z = [float(word) for word in words]
                    if min([x,y,z]) < p_min:
                        p_min = min([x,y,z])
                    if max([x,y,z]) > p_max:
                        p_max = max([x,y,z])
    scale_box = (p_min,p_max)
    
    
    
    
    #=========================Step 1===============================================
    print("Step 1: determine the value of C_expected")
    C_expected, N_C, N_Frames = pa2.compute_C_expected(calbody_filename, calreading_filename)
    D_list, A_list, C_list, N_D, N_A, N_C, N_Frames = pa2.read_calreadings(calreading_filename)
    C = np.transpose(cis.vec_list_to_matrix(C_list))
        
    #=========================Step 2===============================================
    print("Step 2: Produce a suitable distortion correction function")
    #Degree of Berstein Polynomial
    N = 5
    #calculate the bernstein polynomial to fit data
    coefficient = pa2.bernstein_polynomial(C_expected, C, N, scale_box = scale_box )
    # use the correction function
    C_corrected = pa2.correction_function(C.T, coefficient.T, N, scale_box = scale_box )
    
    #=======================Step 3================================================
    print("Step 3: Use the distortion correction function to repeat pivot calibration")
    p_tip_EM, p_pivot_EM, calib_local_frame_EM = pivot_calibration.EM_Pivot_Calibration_With_Correction(empivot_filename, coefficient.T, N, scale_box = scale_box)
    p_tip_EM_d, p_pivot_d_EM, calib_local_frame_EM_d = pivot_calibration.EM_Pivot_Calibration(empivot_filename)
    print('p_pivot without correction\n',p_pivot_d_EM, '\n')
    print('p_pivot after correction\n',p_pivot_EM, '\n')
    
    
    #======================Step 3.1===============================================
    print("Step 3.1: Optical Probe pivot calibration")
    p_tip_OP, p_pivot_OP, calib_local_frame_OP = pivot_calibration.OP_Pivot_Calibration(optpivot_filename, calbody_filename)
    print('p_pivot_Optical:\n',p_tip_OP,sep='')
    
    #=====================Step 4==================================================
    print("Step 4: compute the locations of the fiducials point\n")
    fiducials_2_Base = pa2.fiducials_relative_base(em_fiducials_filename, calib_local_frame_EM, p_tip_EM, coefficient.T, N, scale_box)
    fiducials_2_Base_dist = pa2.fiducials_relative_base(em_fiducials_filename, calib_local_frame_EM_d, p_tip_EM_d, coefficient.T, N, scale_box, correction=False)
    fiducials_2_CT = pa2.fiducials_relative_CT(ct_fiducials_filename)
    
    #=====================Step 5==================================================
    print("Step 5 Compute the registration frame\n")
    F_reg = registration.registration(fiducials_2_Base, fiducials_2_CT)
    F_reg_dist = registration.registration(fiducials_2_Base_dist, fiducials_2_CT)
    # #DEBUG PLOTTER
    # plotter.plot_vec_list(fiducials_2_Base)
    # plotter.plot_vec_list(fiducials_2_CT, color = 'r')
    # for i in range(len(fiducials_2_Base)):
    #     plotter.plot_arrow(fiducials_2_Base[i],fiducials_2_CT[i]-fiducials_2_Base[i])
         
    
    #=====================Step 6==================================================
    print("Step 6 Compute the pointer tip coordinates with respect to the tracker base\n")       
    nav_list = pa2.calc_nav_points(em_nav_filename, calib_local_frame_EM , p_tip_EM, coefficient.T, N, scale_box, F_reg)
    nav_list_dist = pa2.calc_nav_points(em_nav_filename, calib_local_frame_EM_d , p_tip_EM_d, coefficient.T, N, scale_box, F_reg_dist, correction = False )
    
    #Result print out in console
    print("Navigation tool tip without correction:")
    for nav in nav_list_dist:
        print(nav)
    print("Navigation tool tip after correction:")
    for nav in nav_list:
        print(nav)
    
    #Result Comparasion with DEBUG file
    if 'debug' in output2_filename:
        ct_points = []
        with open('../Input Data/{}'.format(output2_filename), 'r') as f:
            for i,line in enumerate(f):
                words = line.split(',')
                #Strip words
                for w in range(len(words)):
                    words[w] = words[w].strip(' .,')
                #Handle file header, containing data info
                if i == 0:
                    continue
                else:
                    #Handle data
                    x, y, z = [float(word) for word in words]
                    p = cis.Vec3D(x,y,z)
                    ct_points.append(p)
        
        
        avg_error_no_correction = sum([ct_points[i].matrix - nav_list_dist[i].matrix for i in range(len(ct_points))])/len(ct_points)
        avg_error = sum([ct_points[i].matrix - nav_list[i].matrix for i in range(len(ct_points)) ])/len(ct_points)
                    
        print("Avg error before correction:\n", avg_error_no_correction, sep='')
        print("Avg error after correction:\n", avg_error, sep='')
    
    #====================Export Output Text File==================================
    
    with open('../OUTPUT/{}'.format(output1_filename), 'w') as f:
        f.write(str(N_C) +', ' + str(N_Frames) + ', ' + output1_filename)
        f.write('\n')
        f.write('  ' + str(round(p_pivot_EM.x, 2)) +',   '+str(round(p_pivot_EM.y, 2)) +',   '+ str(round(p_pivot_EM.z, 2)))
        f.write('\n')
        f.write('  ' + str(round(p_pivot_OP.x,2)) + ',   '+str(round(p_pivot_OP.y, 2)) + ',   ' + str(round(p_pivot_OP.z, 2)))
        f.write('\n')
        for row in range(len(C_expected)):
            f.write('  ' + str(round(C_expected[row, 0], 2)) +',   ' + str(round(C_expected[row, 1], 2)) + ',   ' + str(round(C_expected[row, 2], 2)))
            f.write('\n')
    
    with open('../OUTPUT/{}'.format(output2_filename), 'w') as f:
        f.write(str(len(nav_list)) + ', ' + output2_filename+'\n')
        for nav in nav_list:
            f.write('  ' + str(round(nav.x, 2)) +',   '+str(round(nav.y, 2)) +',   '+ str(round(nav.z, 2)))
            f.write('\n')
                
            


#*****************************************************************************
#Run all input files
#***************************************************************************
if __name__ == '__main__':
    for dataset in 'abcdefghij':
        main(dataset)
            
            
            
            
            
            
            
            
            
            