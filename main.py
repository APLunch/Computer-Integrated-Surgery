#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 22:22:02 2021

@author: ychen215, hfan15
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

#=========================Registration, Problem 4==============================

print("Problem 4: Compute Expected C")

#name = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
name = 'a'

calbody = 'pa1-debug-'+ name +'-calbody.txt'
calreading = 'pa1-debug-'+ name +'-calreadings.txt'
output = 'pa1-debug-'+ name +'-output1.txt'


C_expected = compute_C_expected(calbody, calreading)



average_error, RMSE = compare_c(calbody, calreading, output)




#====================Pivot Calibration, Problem 5 and Problem 6===============

print("Problem 5: EM Problem pivot Calibration")
tp_EM, p_pivot_EM = pivot_calibration.EM_Pivot_Calibration('pa1-debug-a-empivot.txt')
print("tp = \n",tp_EM,'\n\np_pivot=\n',p_pivot_EM,'\n\n',sep='')

print("Problem 6: Optical Probe pivot Calibration")
tp_OP, p_pivot_OP = pivot_calibration.OP_Pivot_Calibration('pa1-debug-a-optpivot.txt', "pa1-debug-a-calbody.txt")
print("tp = \n",tp_OP,'\n\np_pivot=\n',p_pivot_OP,'\n\n',sep='')
