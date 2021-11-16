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
    
    'vector divition'
    def __div__(self, other):
        return Vec3D(self.matrix/other)
    
    'vector right multiplication'
    def __rmul__(self, other):
        return Vec3D(other*self.matrix)
    
    'vector dot product'
    def dot(self,other):
        return np.matmul(self.matrix.T, other.matrix)[0][0]
    
    'get vector norm'
    def norm(self):
        return np.linalg.norm(self.matrix)
    
    
class Rot3D:
    '3D Rotation'
    def __init__(self, matrix):
        if matrix.shape == (3,3):
            self.matrix = matrix
            if abs(np.linalg.det(matrix) - 1) > 0.0001:
                raise Exception("Rot3D Constructor: A Rot3D matrix must have det of 1")
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

def matrix_to_vec_list(mat):
    '''
    Takes in a matrix consists of colonm vectors and return list of corresponding vectors.

    Parameters
    ----------
    mat : numpy array. the matrix consists of colomn vectors

    Returns
    -------
    list(Vec3D)

    '''
    lst=[]
    for vec in mat.T:
        v3d = Vec3D(vec[0],vec[1],vec[2])
        lst.append(v3d)
    
    return lst
    
    
    
    
class Triangle:
    'Triangle class used for mesh'
    def __init__(self, v1, v2, v3):
        self.vertices = (v1,v2,v3)
    
class Mesh:
    'Mesh class'
    def __init__(self):
        self.triangles = dict()
        self.size = 0
    
    
    def add(self,triangle):
        '''
        Add a triangle into the mesh

        Parameters
        ----------
        triangle : cismath.Triangle
            triangle object to be added

        Returns
        -------
        None.

        '''
        self.triangles[self.size] = triangle
        self.size += 1
    
    def bf_closet_pt_on_mesh(self,a):
        '''
        Find the closes point on the mesh triangles, given a point in space 'a'.
        Using brutal force iterative search
        
        Parameters
        ----------
        a : Vec3D
            A point in 3D space

        Returns
        -------
        Vec3D
            The closest point on mesh
        '''
        
        closest_pt = Vec3D(9999999,9999999,9999999)
        #General case
        for (i,tri) in self.triangles.items():
            p,r,q = tri.vertices
            #Configure lstsq
            A = np.hstack([(q-p).matrix, (r-p).matrix])
            b = (a-p).matrix
            res = np.linalg.lstsq(A,b,rcond=None)
            #unpack lambda and mu
            lam = res[0][0][0]
            miu = res[0][1][0]
            c = p+lam*(q-p)+miu*(r-p)
            #Check condition
            if lam >= 0 and miu >= 0 and (lam+miu) <= 1:
                #Closest pt inside triangle
                if (c-a).norm() < (closest_pt-a).norm():
                    closest_pt = c
            else:
                #closest pt on edge
                if lam < 0:
                    c = ProjectOnSegment(c,r,p)
                elif miu < 0:
                    c = ProjectOnSegment(c,p,q)
                elif lam+miu > 1:
                    c = ProjectOnSegment(c,q,r)
                #Replace closest_pt with closer c, if found
                if (c-a).norm() < (closest_pt-a).norm():
                    closest_pt = c
        return closest_pt
                
def ProjectOnSegment(c,p,q):
    '''
    helper function that determines the projection of closest point on a triangle's edge 
        using given parameters

    Returns
    -------
    Vec3D 
        The projected closest point on edge

    '''
    lam = (c-p).dot((q-p))/((q-p).dot(q-p))
    lam_seg = max(0, min(lam,1))
    c_projection = p + lam_seg*(q-p)
    return c_projection
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    