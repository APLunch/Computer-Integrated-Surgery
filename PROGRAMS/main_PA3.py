#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 22:22:02 2021

@author: Yuxin Chen, hfan15
"""
import time
import assignment3_utilities as pa3
import cismath as cis 
from registration import registration
# =============================================================================================
# =======================================Configurations========================================
# =============================================================================================

name = "J"


#first_reading_file = 'PA3-' + name +'-Debug-SampleReadingsTest.txt'

<<<<<<< Updated upstream
first_reading_file = 'PA3-' + name +'-Unknown-SampleReadingsTest.txt'

<<<<<<< Updated upstream
=======
=======
>>>>>>> Stashed changes
name = "F"
first_reading_file = 'PA3-' + name +'-Debug-SampleReadingsTest.txt'
>>>>>>> Stashed changes
output_filename = 'PA3-' + name + '-Output.txt'
debug_output_filename = 'PA3-' + name + '-Debug-Output.txt' 

first_body_filename = 'Problem3-BodyA.txt'

N_A, a_list, a_tip =pa3.read_body(first_body_filename) 

second_body_filename = 'Problem3-BodyB.txt'

N_B, b_list, b_tip =pa3.read_body(second_body_filename)

A_list, B_list, D_list, N_D, N_samples = pa3.read_sample_readings(first_reading_file, N_A, N_B)

d_list = []
for i in range(N_samples):
    A_sublist = A_list[i*N_A:(i+1)*N_A]
    B_sublist = B_list[i*N_B:(i+1)*N_B]
    Fa = registration(a_list, A_sublist)
    Fb = registration(b_list, B_sublist)
    d = Fb.inv()*(Fa * a_tip)
    d_list.append(d)

#Load Mesh
Mesh = pa3.load_mesh_from_file('Problem3Mesh.sur')



# =============================================================================================
# ===================================Brutal Force Approach=====================================
# =============================================================================================
# find c without KDTree
start = time.time() 
c_list = []
for d in d_list:
    #print(Mesh.bf_closet_pt_on_mesh(d))
    c_list.append(Mesh.closest_pt_on_mesh(d))
end = time.time()
print('='*5+'Time Elapsed Using Brutal Force:{:.3f}'.format(end-start)+'='*5) 
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
print(round(average_error_d.x, 3), round(average_error_d.y, 3), round(average_error_d.z, 3))
print('error in c = ')
print(round(average_error_c.x, 3), round(average_error_c.y, 3), round(average_error_c.z, 3))
print('error in e = ')
<<<<<<< Updated upstream
<<<<<<< Updated upstream
print(average_error_e)
print('\n')
"""


# =============================================================================================
# ======================================KD Tree Approach=======================================
# =============================================================================================
# find c with KDTree
start = time.time() 
Mesh.make_tree()
for d in d_list:
    #print(Mesh.bf_closet_pt_on_mesh(d))
    c_list.append(Mesh.closest_pt_on_mesh(d))
end = time.time()
print('='*5+'Time Elapsed Using KDTree:{:.3f}'.format(end-start)+'='*5) 
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
=======
print(round(average_error_e, 3))
>>>>>>> Stashed changes
=======
print(round(average_error_e, 3))
>>>>>>> Stashed changes
    


# =============================================================================================
# =======================================File Output===========================================
# =============================================================================================
    
with open('../OUTPUT/{}'.format(output_filename), 'w') as f:
    f.write(str(N_samples) +' ' + output_filename)
    f.write('\n')
    for row in range(N_samples):
        f.write('  ' + str(round(d_list[row].x, 2)) +',   '+str(round(d_list[row].y, 2)) +',   ' + str(round(d_list[row].z, 2)) +
                ',   ' + str(round(c_list[row].x, 2)) +',   ' + str(round(c_list[row].y, 2)) +',   ' + str(round(c_list[row].z, 2)) +
                ',   ' + str(round(abs((d_list[row] - c_list[row]).norm()), 3)) + ' ')
        f.write('\n')







