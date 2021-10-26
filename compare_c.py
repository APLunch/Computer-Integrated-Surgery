#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 01:47:18 2021

@author: ychen215
"""

import numpy as np
from compute_C_expected import compute_C_expected
from get_C_output import get_C_output

def compare_c(calbody, calreading, output):
    
    C_expected = compute_C_expected(calbody, calreading)

    #get corresponding c_expected from output file
    C_output = get_C_output(output)

    
    #C_expected = np.round(C_expected, 2)

    #b = np.isclose(C_expected, C_output)
    
    RMSE_x = np.sqrt(np.mean((C_expected[:, 0]-C_output[:, 0])**2))
    RMSE_y = np.sqrt(np.mean((C_expected[:, 1]-C_output[:, 1])**2))
    RMSE_z = np.sqrt(np.mean((C_expected[:, 2]-C_output[:, 2])**2))
    
    RMSE = np.sqrt(np.mean((C_expected-C_output)**2))
    
    sum_error_x = 0
    sum_error_y = 0
    sum_error_z = 0
    sum_errow = 0
    
    for i in range(len(C_expected)):
        sum_error_x += abs(C_expected[i, 0]-C_output[i, 0])
        sum_error_y += abs(C_expected[i, 1]-C_output[i, 1])                      
        sum_error_z += abs(C_expected[i, 2]-C_output[i, 2])
    
    sum_error =  sum_error_x + sum_error_y + sum_error_z
    
    average_error_x = sum_error_x/len(C_expected)
    average_error_y = sum_error_y/len(C_expected)
    average_error_z = sum_error_z/len(C_expected)
    average_error = sum_error/len(C_expected)
    
    print("RMSE = ", RMSE)
    print("Average error = ", average_error)
    print('\n')
    
    return average_error, RMSE