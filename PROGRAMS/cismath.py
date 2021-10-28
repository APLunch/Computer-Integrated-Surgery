# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 16:21:15 2021

@author: MikeH
"""

import numpy as np

class Vec3D:
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
            self.x = x
            self.y = y
            self.z = z
            
        elif len(args) == 1:
            matrix = args[0]
            if matrix.shape != (3,1):
                raise Exception("A Vec3D needs to be discribe in a 3x1 matrix(vector)")
            self.matrix = matrix
            self.x = matrix[0,0]
            self.y = matrix[1,0]
            self.z = matrix[2,0]
        
        else:
            raise Exception("Point3D constructor: Wrong number of parameters, 1 or 3 needed, {} given".format(len(args)))
            
    'vector-vector addition'
    def __add__(self, other):
        return Vec3D(self.matrix + other.matrix)

    'vector-vector subtraction'
    def __sub__(self, other):
        return Vec3D(self.matrix - other.matrix)
    
    'vector-scaler multiplication'
    def __mul__(self, other):
        return Vec3D(other*self.matrix)
    
    'vector printout'
    def __str__(self):
        return str(self.matrix)
            
        
    
    
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
        if isinstance(other, Vec3D):
            return Vec3D(np.matmul(self.matrix, other.matrix))
        
        'Matrix - matrix multiplication'
        if isinstance(other, Rot3D):
            return Rot3D(np.matmul(self.matrix, other.matrix))
        
    'Matrix printout'
    def __str__(self):
        return str(self.matrix)
    
    'Inverse matrix'
    def inv(self):
        'return the inversed Rot3D object'
        return Rot3D(np.linalg.inv(self.matrix))

class Frame:
    '3d Frame'
    def __init__(self, R, p):
        self.R = R
        self.p = p
    
    'Frame multiplication'
    def __mul__(self, other):
        if isinstance(other,Frame):
            return Frame(self.R * other.R, self.R * other.p + self.p)
        elif isinstance(other, Vec3D):
            return self.R * other + self.p
        else:
            raise Exception("Frame multiplication type error")
            
    def __str__(self):
        return "====\nR:\n{}\np:\n{}\n====".format(str(self.R),str(self.p))
    
    'Get inversed transformation'
    def inv(self):
        'Returns the inversed rigid transformation'
        return Frame(self.R.inv(),self.R.inv()*(self.p*-1))
    
    
def vec_list_to_matrix(veclist):
    'Takes a list of Vec3D, return the corresponding matrix that contains all colonm vectors'
    m = np.hstack([v.matrix for v in veclist])
    return m


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    