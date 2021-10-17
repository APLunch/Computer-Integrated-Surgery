# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 16:21:15 2021

@author: MikeH
"""

import numpy as np

class Point3D:
    '3D point in space'
    
    def __init__(self, *args):
        if len(args) == 3:
            x = args[0]
            y = args[1]
            z = args[2]
            self.matrix = (np.array(
                                [[x],
                                 [y],
                                 [z]]
                                 ))
        elif len(args) == 1:
            matrix = args[0]
            if matrix.shape != (3,1):
                raise Exception("A Point3D needs to be discribe in a 3x1 matrix(vector)")
            self.matrix = matrix
        
        else:
            raise Exception("Point3D constructor: Wrong number of parameters, 1 or 3 needed, {} given".format(len(args)))
            
    def __add__(self, other):
        return Point3D(self.matrix + other.matrix)
            
        
    
    
class Rot3D:
    '3D Rotation'
    
    def __init__(self, matrix):
        if matrix.shape == (3,3):
            self.matrix = matrix
        else:
            raise Exception("Rot3D Constructor: A Rot3D needs to be discribe in a 3x3 matrix")
            
    def __mul__(self, other):
        'Overload multiplication'
        
        'Matrix-vector multiplication'
        if isinstance(other, Point3D):
            return Point3D(np.matmul(self.matrix, other.matrix))
        
        'Matrix - matrix multiplication'
        if isinstance(other, Rot3D):
            return Rot3D(np.matmul(self.matrix, other.matrix))
        

class Frame:
    '3d Frame'
    def __init__(self, R, p):
        self.R = R
        self.p = p
    
    def __mul__(self, other):
        if isinstance(Frame):
            return Frame(self.R * other.R, self.R * other.p + self.p)
        else:
            raise Exception("Frame multiplication type error")