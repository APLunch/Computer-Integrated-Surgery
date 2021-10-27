#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 01:34:25 2021

@author: Yuxin Chen
"""

from read_calbody import read_calbody
from read_calreadings import read_calreadings
from registration import registration
import numpy as np
import cismath as cis 

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