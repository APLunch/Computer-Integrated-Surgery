# -*- coding: utf-8 -*-
"""
Created on Sun Dec 12 02:21:51 2021

@author: MikeH
"""

import numpy as np
import time
import assignment345_utilities as pa5
import cismath as cis 
from registration import registration
import math
import random

def main(name, bf = False):
    if name in 'ABCDEF':
        reading_file = 'PA5-' + name +'-Debug-SampleReadingsTest.txt'
    else:
        reading_file = 'PA5-' + name +'-Unknown-SampleReadingsTest.txt'
    output_filename = 'PA5-' + name + '-Output.txt'
    debug_output_filename = 'PA5-' + name + '-Debug-Output.txt'
    modes_filename = 'Problem5Modes.txt'
    A_body_filename = 'Problem5-BodyA.txt'
    B_body_filename = 'Problem5-BodyB.txt'
       
    #Load Mesh
    Mesh = pa5.load_mesh_from_file('Problem5MeshFile.sur')
    if not(bf):
        Mesh.make_tree(depth = -1)
        
    #read data files
    N_A, a_list, a_tip =pa5.read_body(A_body_filename)     
    N_B, b_list, b_tip =pa5.read_body(B_body_filename)
    A_list, B_list, D_list, N_D, N_samples = pa5.read_sample_readings(reading_file, N_A, N_B)
    Modes = pa5.read_modes(modes_filename)
        
         
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
        
    # Deformable registration
    F_reg,weights,s_list,c_list, e_list = pa5.ICP_And_Calc_Deformation_Weights(d_list,Mesh,Modes)
    
    
    #read debug output file provided by professor and compare
    if name in 'ABCDEF':
        weights_output, s_output, c_output, e_output = pa5.get_debug_output(debug_output_filename)
        average_error_weights, average_error_s, average_error_c, average_error_e = pa5.compare_output(weights, s_list, c_list, e_list, weights_output, s_output, c_output, e_output, N_samples)
        print('error in weights = ')
        print(average_error_weights)
        print('error in s = ')
        print(round(average_error_s.x, 3), round(average_error_s.y, 3), round(average_error_s.z, 3))
        print('error in c = ')
        print(round(average_error_c.x, 3), round(average_error_c.y, 3), round(average_error_c.z, 3))
        print('error in norm = ')
        print(average_error_e)
        print('\n')
        
    #average_norm = np.mean(e_list)
    #print('average_norm =')
    #print(average_norm)
    
    # =============================================================================================
    # =======================================File Output===========================================
    # =============================================================================================
   
    pa5.output_files_PA5(name, weights, s_list, c_list)
    
    
#Run with all datasets#
for name in 'ABCDEFGHJK':
    print('+++++++++++++++Data Set {}+++++++++++++++++'.format(name))
    main(name)