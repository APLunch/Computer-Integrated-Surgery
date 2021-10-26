#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 22:34:13 2021

@author: ychen215
"""

import numpy as np
import cismath as cis

def read_calreadings(filename):

    #Read files and load EM Calibration data data
    N_D = 0
    N_A = 0
    N_C = 0
    N_Frames = 0
    
    with open("Data/{}".format(filename),'r') as f:
        lines = f.readlines()
        
        #read first line to get ND, NA and NC
        first_line = lines[0].split(",")
        
        N_D = int(first_line[0])
        N_A = int(first_line[1])
        N_C = int(first_line[2])
        N_Frames = int(first_line[3])
        
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