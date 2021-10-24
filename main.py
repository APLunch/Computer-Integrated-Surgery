#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 22:22:02 2021

@author: ychen215
"""

from registration import registration
from read_text import read_text
import numpy as np

d, a, c = read_text('pa1-debug-a-calbody.txt')

D, A, C = read_text('pa1-debug-a-calreadings.txt')


Rd, pd = registration(d, D)
Ra, pa = registration(a, A)
Rc, pc = registration(c, C)
 

C_expected = np.zeros([len(c), 3])
for i in range(len(c)):

    c_vector = np.zeros([1, 3])
    c_vector += c[i]
    
    z = np.matmul(Ra, np.transpose(c_vector)) + pa
    
    C_e = np.matmul(np.transpose(Rd), z) - np.matmul(np.transpose(Rd), pd)
    
    C_expected[i] = np.transpose(C_e)


    
    