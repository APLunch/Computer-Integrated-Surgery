# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 12:59:30 2021

@author: Yunxin Chen
"""

import numpy as np
import cismath as cis 
from registration import registration

# This function read the C value from the given output txt file
def get_C_output(file_name):
    #read txt file line by line
    lines = []
    with open("../Input Data/{}".format(file_name),'r') as f:
        lines = f.readlines()
    f.close()
    
    #read first line to get N_C and N_frames
    first_line = lines[0].split(",")
    N_C = int(first_line[0])
    N_frames = int(first_line[1])
    
    #save the C value from given output txt file as an numpy.ndarray
    C_output = np.zeros([N_C*N_frames, 3])
    for row in range(3, N_C*N_frames+3):
        p = lines[row].split(",")
        C_output[row-3, 0] = float(p[0])
        C_output[row-3, 1] = float(p[1])
        C_output[row-3, 2] = float(p[2])
    
    return C_output


# This function compute the C_expected value as required in Question4
def compute_C_expected(filename1, filename2):
    
    # read the coordinates of a and A 
    d_list, a_list, c_list = read_calbody(filename1)
    D_list, A_list, C_list, N_D, N_A, N_C, N_Frames = read_calreadings(filename2)
    
    # initialize the error in Fa and Fd for unit test
    error_Fa_x = 0
    error_Fa_y = 0
    error_Fa_z = 0
    error_Fd_x = 0
    error_Fd_y = 0
    error_Fd_z = 0

    
    #calculate the C_expected
    C_expected_list = []
    for i in range(N_Frames):
        # calculate the corresponding Fd, Fa, and Fc in each frame
        D_sublist = D_list[i*N_D:(i+1)*N_D]
        A_sublist = A_list[i*N_A:(i+1)*N_A]
        C_sublist = C_list[i*N_C:(i+1)*N_C]
        Fd = registration(d_list, D_sublist)
        Fa = registration(a_list, A_sublist)
        Fc = registration(c_list, C_sublist)
        
    
        # calculate the C_expected by using the corresponding Fd, Fa
        for c_vector in c_list:
            c = Fd.inv()*(Fa * c_vector)
            C_expected_list.append(c)
    
        # Unit test for Fa and Fd by calculating the error in Fa and Fd
        A_test = []
        D_test = []
    
        # calcule late error in Fa
        for n in range(len(a_list)):
            a_vector = a_list[n]
            A_test.append(Fa * a_vector)
            error_Fa_x += abs(A_sublist[n].x - (Fa * a_vector).x)
            error_Fa_y += abs(A_sublist[n].y - (Fa * a_vector).y)
            error_Fa_z += abs(A_sublist[n].z - (Fa * a_vector).z)
            
            
        # calcule late error in Fd
        for n in range(len(d_list)):
            d_vector = d_list[n]
            D_test.append(Fd * d_vector)
            error_Fd_x += abs(D_sublist[n].x - (Fd * d_vector).x)
            error_Fd_y += abs(D_sublist[n].y - (Fd * d_vector).y)
            error_Fd_z += abs(D_sublist[n].z - (Fd * d_vector).z)
     
            
    error_Fa_x = error_Fa_x/len(A_list)
    error_Fa_y = error_Fa_y/len(A_list)
    error_Fa_z = error_Fa_z/len(A_list)
    #print('error Fa = ', error_Fa_x, error_Fa_y, error_Fa_z)
    
    
    error_Fd_x = error_Fd_x/len(D_list)
    error_Fd_y = error_Fd_y/len(D_list)
    error_Fd_z = error_Fd_z/len(D_list)
    #print('error Fd = ', error_Fd_x, error_Fd_y, error_Fd_z)
    print('\n')    
    
    
    print("C_expected:")    
    # Save C_expected as an numpy.ndarray
    C_expected = np.transpose(cis.vec_list_to_matrix(C_expected_list))
    
    print(C_expected)
    print('\n')
    
    return C_expected, N_C, N_Frames


# This function comare the C_expected we calculated and C value from given output file
def compare_c(C_expected, C_output):
    sum_error_x = 0
    sum_error_y = 0
    sum_error_z = 0
    
    for i in range(len(C_expected)):
        sum_error_x += abs(C_expected[i, 0]-C_output[i, 0])
        sum_error_y += abs(C_expected[i, 1]-C_output[i, 1])                      
        sum_error_z += abs(C_expected[i, 2]-C_output[i, 2])
    
    
    #calculate the average error of x, y, z
    average_error_x = sum_error_x/len(C_expected)
    average_error_y = sum_error_y/len(C_expected)
    average_error_z = sum_error_z/len(C_expected)
    
    
    return average_error_x, average_error_y, average_error_z

# read the calbody txt file
def read_body(filename):
    
    with open("../2021 PA 3-5 Student Data/{}".format(filename),'r') as f:
        lines = f.readlines()
        
        #read first line to get ND, NA and NC
        first_line = lines[0].split(" ")
        N_markers = int(first_line[0])
        
        #save all the coordinates of d, a, c to three defined lists of Vec3D object
        points_LED = []
        
        for line in lines[1:N_markers+1]:
            words = line.split()
            words = [word.strip(' ') for word in words]
            print(float(words[0]), float(words[1]), float(words[2]))
            points_LED.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
        
        words = lines[N_markers+1].split()
        words = [word.strip(' ') for word in words]
        print(float(words[0]), float(words[1]), float(words[2]))
        p_tip = cis.Vec3D(float(words[0]), float(words[1]), float(words[2]))
    
    return N_markers, points_LED, p_tip


# read the calreadings txt file
def read_sample_readings(filename, N_A, N_B):
    
    with open("../2021 PA 3-5 Student Data/{}".format(filename),'r') as f:
        lines = f.readlines()
        
        #read first line to get ND, NA and NC
        first_line = lines[0].split(",")
        
        N_S = int(first_line[0])
        N_samples = int(first_line[1])
        
        N_D = N_S - N_A - N_B
        
        #save all the coordinates of A, B, D to three defined lists of Vec3D object

        A_list = []
        B_list = []
        D_list = []
        
        
        for line in lines[1:N_D+1]:
            words = line.split()
            words = [word.strip(',') for word in words]
            points_d.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
    
        for line in lines[1+N_D:N_D+1+N_A]:
            words = line.split()
            words = [word.strip(',') for word in words]
            points_a.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
    
        for line in lines[N_D+1+N_A:N_D+1+N_A+N_C]:
            words = line.split()
            words = [word.strip(',') for word in words]
            points_c.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
        
        """
        for i in range(N_Frames):
            for line in lines[1+i*(N_D+N_A+N_C) : N_D+1+i*(N_D+N_A+N_C)]:
                words = line.split()
                words = [word.strip(',') for word in words]
                D_list.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
        
            for line in lines[1+N_D+i*(N_D+N_A+N_C) : N_D+1+N_A+i*(N_D+N_A+N_C)]:
                words = line.split()
                words = [word.strip(',') for word in words]
                A_list.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
        
            for line in lines[N_D+1+N_A+i*(N_D+N_A+N_C) : N_D+1+N_A+N_C+i*(N_D+N_A+N_C)]:
                words = line.split()
                words = [word.strip(',') for word in words]
                C_list.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
        """
    return D_list, A_list, C_list, N_D, N_A, N_C, N_Frames



def get_calibration_output(filename):
    with open("../Input Data/{}".format(filename),'r') as f:
        lines = f.readlines()
        #Get output tip position for em probe
        pt_em_x, pt_em_y, pt_em_z = [float(word.strip(' ,.')) for word in lines[1].split()]
        pt_em = cis.Vec3D(pt_em_x, pt_em_y, pt_em_z)
        
        #Get output tip position for optical probe
        pt_op_x, pt_op_y, pt_op_z = [float(word.strip(' ,.')) for word in lines[2].split()]
        pt_op = cis.Vec3D(pt_op_x, pt_op_y, pt_op_z)
        
        return pt_em, pt_op
        
        
        
        
        
        
        
        
        
        
        
        
        
        