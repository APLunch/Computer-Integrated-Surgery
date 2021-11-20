#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 22:22:02 2021

@author: Yuxin Chen, hfan15
"""
import numpy as np
import time
import assignment3_utilities as pa3
import cismath as cis 
from registration import registration

def main(name):
    # =============================================================================================
    # =======================================Configurations========================================
    # =============================================================================================
    if name in 'ABCDEF':
        reading_file = 'PA3-' + name +'-Debug-SampleReadingsTest.txt'
    else:
        reading_file = 'PA3-' + name +'-Unknown-SampleReadingsTest.txt'
    
    
    output_filename = 'PA3-' + name + '-Output.txt'
    
    A_body_filename = 'Problem3-BodyA.txt'
    B_body_filename = 'Problem3-BodyB.txt'
    
    
    #read data files
    N_A, a_list, a_tip =pa3.read_body(A_body_filename) 
    
    N_B, b_list, b_tip =pa3.read_body(B_body_filename)
    
    A_list, B_list, D_list, N_D, N_samples = pa3.read_sample_readings(reading_file, N_A, N_B)
    
    
    # initialize F0
    F0 = cis.Frame(cis.Rot3D(np.eye(3)),cis.Vec3D(0,0,0))
    
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
    
    #Load Mesh
    Mesh = pa3.load_mesh_from_file('Problem3MeshFile.sur')
    
    
    """
    # =============================================================================================
    # ===================================Brutal Force Approach=====================================
    # =============================================================================================
    # find c without KDTree
    start = time.time() 
    c_list = []
    for d in d_list:
        #print(Mesh.bf_closet_pt_on_mesh(d))
        c_list.append(Mesh.closest_pt_on_mesh(F0*d))
    end = time.time()
    print('='*5+'Time Elapsed Using Brutal Force: {:.3f}s'.format(end-start)+'='*5) 
    # calculate norm of d minus c
    e = []
    for row in range(N_samples):
        e.append(abs((d_list[row] - c_list[row]).norm()))
    
    
    
    #read debug output file provided by professor
    d_output, c_output, e_output = pa3.get_debug_output(debug_output_filename)
    # compare our output with debug output
    average_error_d, average_error_c, average_error_e = pa3.compare_output(d_list, c_list, e, d_output, c_output, e_output, N_samples)
    print('error in d = ')
    print(round(average_error_d.x, 3), round(average_error_d.y, 3), round(average_error_d.z, 3))
    print('error in c = ')
    print(round(average_error_c.x, 3), round(average_error_c.y, 3), round(average_error_c.z, 3))
    print('error in e = ')
    print(average_error_e)
    print('\n')
    """
    
    
    # =============================================================================================
    # ======================================KD Tree Approach=======================================
    # =============================================================================================
    # find c with KDTree
    c_list = [] 
    start = time.time() 
    Mesh.make_tree()
    for d in d_list:
        #print(Mesh.bf_closet_pt_on_mesh(d))
        c_list.append(Mesh.closest_pt_on_mesh(F0*d))
    end = time.time()
    print('='*5+'Time Elapsed Using KDTree: {:.3f}s'.format(end-start)+'='*5) 
    # calculate e
    e = []
    for row in range(N_samples):
        e.append(abs((d_list[row] - c_list[row]).norm()))
    
    """
    #read debug output file provided by professor
    d_output, c_output, e_output = pa3.get_debug_output(debug_output_filename)
    # compare our output with debug output
    average_error_d, average_error_c, average_error_e = pa3.compare_output(d_list, c_list, e, d_output, c_output, e_output, N_samples)
    print('error in d = ')
    print(round(average_error_d.x,3), round(average_error_d.y, 3), round(average_error_d.z, 3))
    print('error in c = ')
    print(round(average_error_c.x, 3), round(average_error_c.y, 3), round(average_error_c.z, 3))
    print('error in e = ')
    print(round(average_error_e, 3))
    
    """
    
    
    # =============================================================================================
    # =======================================File Output===========================================
    # =============================================================================================
        
    with open('../OUTPUT/{}'.format(output_filename), 'w') as f:
        f.write(str(N_samples) +' ' + output_filename)
        f.write('\n')
        for row in range(N_samples):
            s = '  {:9.2f}{:9.2f}{:9.2f}   {:9.2f}{:9.2f}{:9.2f}   {:9.3f}'.format(d_list[row].x,d_list[row].y,d_list[row].z,c_list[row].x,
                                                                             c_list[row].y,c_list[row].z,(d_list[row] - c_list[row]).norm())
            f.write(s)
            
            f.write('\n')



#Call main for each input file
if __name__ == '__main__':
    for name in 'ABCDEFGHJ':
        print("++++++++++++++++DataSet {}+++++++++++++++++++++++".format(name))
        main(name)
        print("-------------------------------------------------\n")



