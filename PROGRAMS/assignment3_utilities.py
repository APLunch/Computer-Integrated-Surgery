# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 12:59:30 2021

@author: Yunxin Chen
"""

import numpy as np
import cismath as cis 
from registration import registration



# read the calbody txt file
def read_body(filename):
    
    with open("../2021 PA 3-5 Student Data/{}".format(filename),'r') as f:
        lines = f.readlines()
        
        #read first line to get ND, NA and NC
        first_line = lines[0].split(" ")
        N_markers = int(first_line[0])
        
        #save all the coordinates of d, a, c to three defined lists of Vec3D object
        points_LED = []
        
        for line in lines[1:N_markers+1]:
            words = line.split()
            words = [word.strip(' ') for word in words]
            print(float(words[0]), float(words[1]), float(words[2]))
            points_LED.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
        
        words = lines[N_markers+1].split()
        words = [word.strip(' ') for word in words]
        print(float(words[0]), float(words[1]), float(words[2]))
        p_tip = cis.Vec3D(float(words[0]), float(words[1]), float(words[2]))
    
    return N_markers, points_LED, p_tip


# read the calreadings txt file
def read_sample_readings(filename, N_A, N_B):
    
    with open("../2021 PA 3-5 Student Data/{}".format(filename),'r') as f:
        lines = f.readlines()
        
        #read first line to get ND, NA and NC
        first_line = lines[0].split(",")
        
        N_S = int(first_line[0])
        N_samples = int(first_line[1])
        
        N_D = N_S - N_A - N_B
        
        #save all the coordinates of A, B, D to three defined lists of Vec3D object

        A_list = []
        B_list = []
        D_list = []
        
        for i in range(N_samples):
            for line in lines[1+i*(N_S) : N_A+1+i*(N_S)]:
                words = line.split()
                words = [word.strip(',') for word in words]
                A_list.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
        
                print(float(words[0]), float(words[1]), float(words[2]))
        
            for line in lines[1+N_A+i*(N_S) : N_A+1+N_B+i*(N_S)]:
                words = line.split()
                words = [word.strip(',') for word in words]
                B_list.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
                
                print(float(words[0]), float(words[1]), float(words[2]))
                
            for line in lines[N_A+1+N_B+i*(N_S) : N_A+1+N_B+N_D+i*(N_S)]:
                words = line.split()
                words = [word.strip(',') for word in words]
                D_list.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
        
                print(float(words[0]), float(words[1]), float(words[2]))
                
    return A_list, B_list, D_list, N_D, N_samples


        
        
        
        
        
        
        
        
        
        
        
        
        
        