#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 00:21:17 2021

@author: ychen215
"""


import numpy as np
import cismath


def registration(input_points_a, input_points_b):
    '''
    This function accepts two lists of Vec3D object, and returns the transformation between the two sets of vectors.
    '''
    points_a = np.zeros([len(input_points_a),3])
    points_b = np.zeros([len(input_points_b),3])
    for row,pt in enumerate(input_points_a):
        points_a[row,0] = pt.x
        points_a[row,1] = pt.y
        points_a[row,2] = pt.z
    for row,pt in enumerate(input_points_b):
        points_b[row,0] = pt.x
        points_b[row,1] = pt.y
        points_b[row,2] = pt.z
    #====================以上是我把传入的参数改成了Vec3D列表，以下依然用以前的格式就行==========================#
    
    a_sum = np.zeros([1, 3])
    b_sum = np.zeros([1, 3])

    for i in range(len(points_a)):
        a_sum += points_a[i]
        b_sum += points_b[i]

    a_mean = a_sum/len(points_a)
    b_mean = b_sum/len(points_b)
    
    a_p = points_a - a_mean
    b_p = points_b - b_mean
    

    M = None

    for i in range(len(points_a)):
        m = np.array([[0.0,                   b_p[i][0]-a_p[i][0],    b_p[i][1]-a_p[i][1],    b_p[i][2]-a_p[i][2]],
                      [b_p[i][0]-a_p[i][0], 0.0,                      -(b_p[i][2]+a_p[i][2]), b_p[i][1]+a_p[i][1]],
                      [b_p[i][1]-a_p[i][1], b_p[i][2]+a_p[i][2],    0.0,                      -(b_p[i][0]+a_p[i][0])],
                      [b_p[i][2]-a_p[i][2], -(b_p[i][1]+a_p[i][1]), b_p[i][0]+a_p[i][0],    0.0]])
        
        if i == 0:
            M = m
        else:
            M = np.append(M, m, axis=0)
    
    [u, s, vh] = np.linalg.svd(M)
    #print(vh)
    q = vh[3, :]
            
    R = np.array([[q[0]**2+q[1]**2-q[2]**2-q[3]**2, 2*(q[1]*q[2]-q[0]*q[3]),         2*(q[1]*q[3]+q[0]*q[2])],
                  [2*(q[1]*q[2]+q[0]*q[3]),         q[0]**2-q[1]**2+q[2]**2-q[3]**2, 2*(q[2]*q[3]-q[0]*q[1])],
                  [2*(q[1]*q[3]-q[0]*q[2]),         2*(q[2]*q[3]+q[0]*q[1]),         q[0]**2-q[1]**2-q[2]**2+q[3]**2]])
        
    #print(np.linalg.det(R))
            
    p = np.transpose(b_mean) - np.matmul(R, np.transpose(a_mean))
    '''  
    print('R=')
    print(R)
    print('P=')
    print(p)
    '''
    
    return cismath.Frame(cismath.Rot3D(R), cismath.Vec3D(p))


"""

H = np.zeros([3, 3])
for i in range(len(points_a)):
    H += [[a_p[i][0]*b_p[i][0], a_p[i][0]*b_p[i][1], a_p[i][0]*b_p[i][2]],
          [a_p[i][1]*b_p[i][0], a_p[i][1]*b_p[i][1], a_p[i][1]*b_p[i][2]],
          [a_p[i][2]*b_p[i][0], a_p[i][2]*b_p[i][1], a_p[i][2]*b_p[i][2]]]


print(H)

[u, s, vh] = np.linalg.svd(H, full_matrices=True)
    
R = np.transpose(vh) * np.transpose(u)

print(R)
    
if np.linalg.det(R) != 1:
    print('error')
    print(np.linalg.det(R))
        
p = b_mean - R * a_mean


"""