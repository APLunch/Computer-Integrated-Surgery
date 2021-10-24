#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 14:48:05 2021

@author: Yuxin(Ethan) Chen
"""


import numpy as np


def read_text(file_name):

    #read txt file line by line
    lines = []
    with open(file_name, "r") as f:
        lines = f.readlines()
    f.close()
    
    #read first line to get ND, NA and NC
    first_line = lines[0].split(",")
    
    N_D = int(first_line[0])
    N_A = int(first_line[1])
    N_C = int(first_line[2])
    
    
    d = np.zeros([N_D, 3])
    for row in range(1, N_D+1):
        p = lines[row].split(",")
        d[row-1, 0] = float(p[0])
        d[row-1, 1] = float(p[1])
        d[row-1, 2] = float(p[2])
        
    
    a = np.zeros([N_A, 3])
    for row in range(N_D+1, N_D+N_A+1):
        p = lines[row].split(",")
        a[row -1 - N_D, 0] = float(p[0])
        a[row -1 - N_D, 1] = float(p[1])
        a[row -1 - N_D, 2] = float(p[2])
        
    
    c = np.zeros([N_C, 3])
    for row in range(N_D+N_A+1, N_D+N_A+N_C+1):
        p = lines[row].split(",")
        c[row -1 - N_D - N_A, 0] = float(p[0])
        c[row -1 - N_D - N_A, 1] = float(p[1])
        c[row -1 - N_D - N_A, 2] = float(p[2])
     
    
    
    return d, a, c