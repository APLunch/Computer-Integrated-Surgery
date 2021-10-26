#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 21:22:02 2021

@author: ychen215
"""
import numpy as np

def get_C_output(file_name):
    #read txt file line by line
    lines = []
    with open("Data/{}".format(file_name),'r') as f:
        lines = f.readlines()
    f.close()
    
    #read first line to get ND, NA and NC
    first_line = lines[0].split(",")
    
    N_C = int(first_line[0])
    N_frames = int(first_line[1])
    
    #C = np.zeros([N_C, 3])    
    C_output = np.zeros([N_C*N_frames, 3])
    for row in range(3, N_C*N_frames+3):
        p = lines[row].split(",")
        C_output[row-3, 0] = float(p[0])
        C_output[row-3, 1] = float(p[1])
        C_output[row-3, 2] = float(p[2])
    
    return C_output