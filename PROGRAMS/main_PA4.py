#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 22:22:02 2021

@author: Yuxin Chen, hfan15
"""
import numpy as np
import time
import assignment3_utilities as pa4
import cismath as cis 
from registration import registration
import math
import random

def main(name, bf = False):  
    if name in 'ABCDEF':
        reading_file = 'PA4-' + name +'-Debug-SampleReadingsTest.txt'
    else:
        reading_file = 'PA4-' + name +'-Unknown-SampleReadingsTest.txt'
    output_filename = 'PA4-' + name + '-Output.txt'
    debug_output_filename = 'PA4-' + name + '-Debug-Output.txt'
    A_body_filename = 'Problem4-BodyA.txt'
    B_body_filename = 'Problem4-BodyB.txt'
       
    #Load Mesh
    Mesh = pa4.load_mesh_from_file('Problem4MeshFile.sur')
    if not(bf):
        Mesh.make_tree(depth = -1)
        
    #read data files
    N_A, a_list, a_tip =pa4.read_body(A_body_filename)     
    N_B, b_list, b_tip =pa4.read_body(B_body_filename)
    A_list, B_list, D_list, N_D, N_samples = pa4.read_sample_readings(reading_file, N_A, N_B)
        
         
    # calculate d
    d_list = []
    for i in range(N_samples):
        A_sublist = A_list[i*N_A:(i+1)*N_A]
        B_sublist = B_list[i*N_B:(i+1)*N_B]
        #calculate Fa
        Fa = registration(a_list, A_sublist)
        #calculate Fb
        Fb = registration(b_list, B_sublist)
            
        #calculate d
        d = Fb.inv()*(Fa * a_tip)
        d_list.append(d)
        
    
    #ICP
    start = time.time()    
    s_list, c_list, e_list = pa4.ICP(d_list, Mesh)
    end = time.time()
    print('='*5+'Time Elapsed Using KDTree: {:.3f}s'.format(end-start)+'='*5) 
    
    
    #read debug output file provided by professor
    if name in 'ABCDEF':
        s_output, c_output, e_output = pa4.get_debug_output(debug_output_filename)
        # compare our output with debug output
        average_error_s, average_error_c, average_error_e = pa4.compare_output(s_list, c_list, e_list, s_output, c_output, e_output, N_samples)
        print('error in s = ')
        print(round(average_error_s.x, 3), round(average_error_s.y, 3), round(average_error_s.z, 3))
        print('error in c = ')
        print(round(average_error_c.x, 3), round(average_error_c.y, 3), round(average_error_c.z, 3))
        print('error in norm = ')
        print(average_error_e)
        print('\n')
        
    
    
    # =============================================================================================
    # =======================================File Output===========================================
    # =============================================================================================
    with open('../OUTPUT/{}'.format(output_filename), 'w') as f:
        f.write(str(N_samples) +' ' + output_filename)
        f.write('\n')
        for row in range(N_samples):
            sp = '  {:9.2f}{:9.2f}{:9.2f}   {:9.2f}{:9.2f}{:9.2f}   {:9.3f}'.format(s_list[row].x, s_list[row].y, s_list[row].z, c_list[row].x,
                 c_list[row].y,c_list[row].z, e_list[row])
            f.write(sp)
                
            f.write('\n')
    icp_time = end-start
    return icp_time


#Run with all datasets#
for name in 'ABCDEFGHJ':
    print('+++++++++++++++Data Set {}+++++++++++++++++'.format(name))
    main(name)



