#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 22:22:02 2021

@author: ychen215, hfan15
"""

from registration import registration
from compute_C_expected import compute_C_expected
from get_C_output import get_C_output
from read_calreadings import read_calreadings
from read_calbody import read_calbody
import numpy as np
import cismath
import pivot_calibration
import cismath as cis

#=========================Registration, Problem 4==============================

C_expected = compute_C_expected('pa1-debug-a-calbody.txt', 'pa1-debug-a-calreadings.txt')

#get corresponding c_expected from output file
C_output = get_C_output("Data/pa1-debug-a-output1.txt")

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



#====================Pivot Calibration, Problem 5 and Problem 6===============

print("Problem 5: EM Problem pivot Calibration")
tp_EM, p_pivot_EM = pivot_calibration.EM_Pivot_Calibration('pa1-debug-a-empivot.txt')
print("tp = \n",tp_EM,'\n\np_pivot=\n',p_pivot_EM,'\n\n',sep='')

print("Problem 6: Optical Probe pivot Calibration")
tp_OP, p_pivot_OP = pivot_calibration.OP_Pivot_Calibration('pa1-debug-a-optpivot.txt', "pa1-debug-a-calbody.txt")
print("tp = \n",tp_OP,'\n\np_pivot=\n',p_pivot_OP,'\n\n',sep='')
