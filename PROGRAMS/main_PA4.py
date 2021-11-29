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
import math
import plotter

name = "B"

reading_file = 'PA4-' + name +'-Debug-SampleReadingsTest.txt'
output_filename = 'PA4-' + name + '-Output.txt'
debug_output_filename = 'PA4-' + name + '-Debug-Output.txt'

A_body_filename = 'Problem4-BodyA.txt'
B_body_filename = 'Problem4-BodyB.txt'
   
#Load Mesh
Mesh = pa3.load_mesh_from_file('Problem4MeshFile.sur')
    
#read data files
N_A, a_list, a_tip =pa3.read_body(A_body_filename)     
N_B, b_list, b_tip =pa3.read_body(B_body_filename)
A_list, B_list, D_list, N_D, N_samples = pa3.read_sample_readings(reading_file, N_A, N_B)
    
     
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
    


# Step 0: Initialization

F = []
# initialize F0
F.append(cis.Frame(cis.Rot3D(np.eye(3)),cis.Vec3D(0,0,0)))

threshold = []
# initialize threshold0
threshold.append(5.0)

sigma = []
sigma.append(0.0)
error_max = []
error_max.append(0.0)
error_mean = []
error_mean.append(1.0)

check = 0

n = 0

start = time.time() 
Mesh.make_tree(depth = -1)
#Input: M, d, F0, n0 check



while check < 10:
    # step 1: matching
    
    # create two empty lists A and B, which is different from Body A and Body B
    A_group = []
    B_group = []
    
    
    
    
    c_list = []
    s_list = []
    e_list = []
    
    for d in d_list:
        #print(Mesh.bf_closet_pt_on_mesh(d))
        c = Mesh.closest_pt_on_mesh(F[n] * d)
        c_list.append(c)
        
        s = F[n]*d
        s_list.append(s)
        
        e = (s -c).norm()
        e_list.append(e)
        
        if e < threshold[n]:
            A_group.append(d)
            B_group.append(c)
     
            
    # step 2: transformation part
    n += 1
    F.append(registration(A_group, B_group))
    
    
    E_list = []
    for k in range(len(A_group)):
        E = B_group[k] - F[n] * A_group[k]
        E_list.append(E)
    
    
    #calculate sigma, error_max and mean error
    sum_dot = 0
    sum_root = 0
    em = 0 
    for k in range(len(E_list)):
        e_dot = E_list[k].dot(E_list[k])
            
        sum_dot += e_dot
        
        if em < math.sqrt(e_dot):
            em = math.sqrt(e_dot)
            
        sum_root += math.sqrt(e_dot)
        
    sigma.append(math.sqrt(sum_dot)/len(E_list))
    error_max.append(em)
    error_mean.append(sum_root/len(E_list))
    
    
    
    # Step 3: adjustment
    
    #adjust threshold
    threshold.append(3*error_mean[n])
    
    
    # Step 4 :(iteration) Termination 
    
    if 0.98 <= error_mean[n]/error_mean[n-1] <= 1 :
        check += 1
    else:
        check = 0
    




end = time.time()
print('='*5+'Time Elapsed Using KDTree: {:.3f}s'.format(end-start)+'='*5) 


 

#read debug output file provided by professor
s_output, c_output, e_output = pa3.get_debug_output(debug_output_filename)
# compare our output with debug output
average_error_s, average_error_c, average_error_e = pa3.compare_output(s_list, c_list, e_list, s_output, c_output, e_output, N_samples)
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




"""
# step 1: matching

# create two empty lists A and B, which is different from Body A and Body B
A_group = []
B_group = []


start = time.time() 
Mesh.make_tree()

c_list = []
e_list = []

for d in d_list:
    #print(Mesh.bf_closet_pt_on_mesh(d))
    c = Mesh.closest_pt_on_mesh(F[n] * d)
    c_list.append(c)
    e = (c - F[n]*d).norm()
    e_list.append(e)
    
    if e < threshold[n]:
        A_group.append(d)
        B_group.append(c)
 
        
# step 2: transformation part

n += 1
F.append(registration(A_group, B_group))

E_list = []
for k in range(len(A_group)):
    E = B_group[k] - F[n] * A_group[k]
    E_list.append(E)


#calculate sigma, error_max and mean error
sum_dot = 0
sum_root = 0
em = 0 
for k in range(len(E_list)):
    e_dot = E_list[k].dot(E_list[k])
        
    sum_dot += e_dot
    
    if em < math.sqrt(e_dot):
        em = math.sqrt(e_dot)
        
    sum_root += math.sqrt(e_dot)
    
sigma.append(math.sqrt(sum_dot)/len(E_list))
error_max.append(em)
error_mean.append(sum_root/len(E_list))



# Step 3: adjustment

#adjust threshold
threshold.append(3*error_mean[n])


if 0.95 <= error_mean[n]/error_mean[n-1] <= 1:
    check += 1
else:
    check = 0


"""

"""
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
    Mesh = pa3.load_mesh_from_file('Problem3Mesh.sur')
    
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
    
    
    
    # =============================================================================================
    # ======================================KD Tree Approach=======================================
    # =============================================================================================
    # find c with KDTree
    start = time.time() 
    c_list = []
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


"""
