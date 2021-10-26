#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 22:22:02 2021

@author: ychen215, hfan15
"""

from registration import registration
from read_text import read_text
import numpy as np
import cismath
import pivot_calibration

#=========================Registration, Problem 4==============================
print("Problem 4: Expected C")
d, a, c = read_text('pa1-debug-a-calbody.txt')

D, A, C = read_text('pa1-debug-a-calreadings.txt')

d_list = []
D_list = []
c_list = []
C_list = []
a_list = []
A_list = []

for row in d:
    d_list.append(cismath.Vec3D(row[0],row[1],row[2]))
for row in D:
    D_list.append(cismath.Vec3D(row[0],row[1],row[2]))
for row in c:
    c_list.append(cismath.Vec3D(row[0],row[1],row[2]))
for row in C:
    C_list.append(cismath.Vec3D(row[0],row[1],row[2]))
for row in a:
    a_list.append(cismath.Vec3D(row[0],row[1],row[2]))
for row in A:
    A_list.append(cismath.Vec3D(row[0],row[1],row[2]))

Fd = registration(d_list, D_list)
Fa = registration(a_list, A_list)
Fc = registration(c_list, C_list)

C_expected = np.zeros([len(c), 3])

print("Expected C:")
for i in range(len(c)):

    c_vector = np.zeros([1, 3])
    c_vector += c[i]
    
    z = np.matmul(Fa.R.matrix, np.transpose(c_vector)) + Fa.p.matrix
    
    C_e = np.matmul(np.transpose(Fd.R.matrix), z) - np.matmul(np.transpose(Fd.R.matrix), Fd.p.matrix)
    
    C_expected[i] = np.transpose(C_e)
    print(C_expected[i])
print('\n')
#====================Pivot Calibration, Problem 5 and Problem 6===============
print("Problem 5: EM Problem pivot Calibration")
tp_EM, p_pivot_EM = pivot_calibration.EM_Pivot_Calibration('pa1-debug-a-empivot.txt')
print("tp = \n",tp_EM,'\n\np_pivot=\n',p_pivot_EM,'\n\n',sep='')

print("Problem 6: Optical Probe pivot Calibration")
tp_OP, p_pivot_OP = pivot_calibration.OP_Pivot_Calibration('pa1-debug-a-optpivot.txt', "pa1-debug-a-calbody.txt")
print("tp = \n",tp_OP,'\n\np_pivot=\n',p_pivot_OP,'\n\n',sep='')