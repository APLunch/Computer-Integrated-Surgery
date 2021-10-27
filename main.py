#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 22:22:02 2021

@author: Yuxin Chen, hfan15
"""

from registration import registration
from compute_C_expected import compute_C_expected
from get_C_output import get_C_output
from compare_c import compare_c
from read_calreadings import read_calreadings
from read_calbody import read_calbody
import numpy as np
import cismath
import pivot_calibration
import cismath as cis


name = 'a'

calbody_filename = 'pa1-debug-'+ name +'-calbody.txt'
calreading_filename = 'pa1-debug-'+ name +'-calreadings.txt'
empivot_filename = 'pa1-debug-' + name + '-empivot.txt'
optpivot_filename = 'pa1-debug-' + name +'-optpivot.txt'
output_filename = 'pa1-debug-'+ name +'-output1.txt'

"""
calbody_filename = 'pa1-unknown-'+ name +'-calbody.txt'
calreading_filename = 'pa1-unknown-'+ name +'-calreadings.txt'
empivot_filename = 'pa1-unknown-' + name + '-empivot.txt'
optpivot_filename = 'pa1-unknown-' + name +'-optpivot.txt'
"""


#=========================Registration, Problem 4==============================

print("Problem 4: Compute Expected C")
C_expected, N_C, N_Frames = compute_C_expected(calbody_filename, calreading_filename)

# read the C values from given output file
C_output = get_C_output(output_filename)
# Compare the C_expected value with C values from given output file
average_error, average_error_x, average_error_y, average_error_z = compare_c(C_expected, C_output)

print('Average error = ')
print(average_error, average_error_x, average_error_y, average_error_z)
print('\n')


#====================Pivot Calibration, Problem 5 and Problem 6===============

print("Problem 5: EM Problem pivot Calibration")
tp_EM, p_pivot_EM = pivot_calibration.EM_Pivot_Calibration(empivot_filename)
print("tp = \n",tp_EM,'\n\np_pivot=\n',p_pivot_EM,'\n\n',sep='')

print("Problem 6: Optical Probe pivot Calibration")
tp_OP, p_pivot_OP = pivot_calibration.OP_Pivot_Calibration(optpivot_filename, calbody_filename)
print("tp = \n",tp_OP,'\n\np_pivot=\n',p_pivot_OP,'\n\n',sep='')


#====================Export Output Text File==================================

with open('pa1-debug-'+ name +'-output1.txt', 'w') as f:
    f.write(str(N_C) +', ' + str(N_Frames) + ', ' + 'pa1-debug-'+ name +'-output1.txt')
    f.write('\n')
    f.write('  ' + str(round(tp_EM.x, 2)) +',   '+str(round(tp_EM.y, 2)) +',   '+ str(round(tp_EM.z, 2)))
    f.write('\n')
    f.write('  ' + str(round(tp_OP.x,2)) + ',   '+str(round(tp_OP.y, 2)) + ',   ' + str(round(tp_OP.z, 2)))
    f.write('\n')

    for row in range(len(C_expected)):
        f.write('  ' + str(round(C_expected[row, 0], 2)) +',   ' + str(round(C_expected[row, 1], 2)) + ',   ' + str(round(C_expected[row, 2], 2)))
        f.write('\n')