#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 23:34:22 2021

@author: ychen215
"""

import cismath as cis

def read_calbody(filename):
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
        #N_frames = int(first_line[3])
        
        points_d = []
        points_a = []
        points_c = []
        
        for line in lines[1:N_D+1]:
            words = line.split()
            words = [word.strip(',') for word in words]
            points_d.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
    
        for line in lines[1+N_D:N_D+1+N_A]:
            words = line.split()
            words = [word.strip(',') for word in words]
            points_a.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
    
        for line in lines[N_D+1+N_A:N_D+1+N_A+N_C]:
            words = line.split()
            words = [word.strip(',') for word in words]
            points_c.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
    
    return points_d, points_a, points_c
    
    
    