# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 12:59:30 2021

@author: Yunxin Chen
"""

import numpy as np
import cismath as cis 
from registration import registration

def get_C_output(file_name):
    #read txt file line by line
    lines = []
    with open("Data/{}".format(file_name),'r') as f:
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

def compute_C_expected(filename1, filename2):
    
    d_list, a_list, c_list = read_calbody(filename1)
    D_list, A_list, C_list, N_D, N_A, N_C, N_Frames = read_calreadings(filename2)
    
    print("C_expected:")
    
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
            
    # Save C_expected as an numpy.ndarray
    C_expected = np.transpose(cis.vec_list_to_matrix(C_expected_list))
    
    print(C_expected)
    print('\n')
    
    return C_expected, N_C, N_Frames



def compare_c(C_expected, C_output):
    
    #C_expected = compute_C_expected(calbody, calreading)

    #get corresponding c_expected from output file
    #C_output = get_C_output(output)

    
    #C_expected = np.round(C_expected, 2)

    #b = np.isclose(C_expected, C_output)
    
    #RMSE_x = np.sqrt(np.mean((C_expected[:, 0]-C_output[:, 0])**2))
    #RMSE_y = np.sqrt(np.mean((C_expected[:, 1]-C_output[:, 1])**2))
    #RMSE_z = np.sqrt(np.mean((C_expected[:, 2]-C_output[:, 2])**2))
    
    #RMSE = np.sqrt(np.mean((C_expected-C_output)**2))
    
    sum_error_x = 0
    sum_error_y = 0
    sum_error_z = 0
    sum_error = 0
    
    for i in range(len(C_expected)):
        sum_error_x += abs(C_expected[i, 0]-C_output[i, 0])
        sum_error_y += abs(C_expected[i, 1]-C_output[i, 1])                      
        sum_error_z += abs(C_expected[i, 2]-C_output[i, 2])
    
    sum_error =  sum_error_x + sum_error_y + sum_error_z
    
    average_error_x = sum_error_x/len(C_expected)
    average_error_y = sum_error_y/len(C_expected)
    average_error_z = sum_error_z/len(C_expected)
    average_error = sum_error/len(C_expected)
    
    #print("RMSE = ", RMSE)
    #print("Average error = ", average_error)
    print('\n')
    
    return average_error_x, average_error_y, average_error_z


def read_calbody(filename):
    
    with open("Data/{}".format(filename),'r') as f:
        lines = f.readlines()
        
        #read first line to get ND, NA and NC
        first_line = lines[0].split(",")
        N_D = int(first_line[0])
        N_A = int(first_line[1])
        N_C = int(first_line[2])
        
        #save all the coordinates of d, a, c to three defined lists of Vec3D object
        points_d = []
        points_a = []
        points_c = []
        
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
    
    return points_d, points_a, points_c


def read_calreadings(filename):
    
    with open("Data/{}".format(filename),'r') as f:
        lines = f.readlines()
        
        #read first line to get ND, NA and NC
        first_line = lines[0].split(",")
        
        N_D = int(first_line[0])
        N_A = int(first_line[1])
        N_C = int(first_line[2])
        N_Frames = int(first_line[3])
        
        #save all the coordinates of D, A, C to three defined lists of Vec3D object
        D_list = []
        A_list = []
        C_list = []
        
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
        
    return D_list, A_list, C_list, N_D, N_A, N_C, N_Frames



def get_calibration_output(filename):
    with open("Data/{}".format(filename),'r') as f:
        lines = f.readlines()
        #Get output tip position for em probe
        pt_em_x, pt_em_y, pt_em_z = [float(word.strip(' ,.')) for word in lines[1].split()]
        pt_em = cis.Vec3D(pt_em_x, pt_em_y, pt_em_z)
        
        #Get output tip position for optical probe
        pt_op_x, pt_op_y, pt_op_z = [float(word.strip(' ,.')) for word in lines[2].split()]
        pt_op = cis.Vec3D(pt_op_x, pt_op_y, pt_op_z)
        
        return pt_em, pt_op
        
        
        
        
        
        
        
        
        
        
        
        
        
        