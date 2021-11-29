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
    def __truediv__(self, other):
        return Vec3D(self.matrix/other)
    
    'vector right multiplication'
    def __rmul__(self, other):
        return Vec3D(other*self.matrix)
    
    'vector dot product'
    def dot(self,other):
        return np.matmul(self.matrix.T, other.matrix)[0][0]
    
    'vector cross product'
    def cross(self, other):
        cp = np.cross(self.matrix.T, other.matrix.T)[0]
        return Vec3D(cp[0], cp[1], cp[2])
    
    'get vector norm'
    def norm(self):
        return np.linalg.norm(self.matrix)
    
    'equality'
    def __eq__(self,other):
        return np.isclose(self.x, other.x) and np.isclose(self.y, other.y) and np.isclose(self.z, other.z)
    
    'hash func'
    def __hash__(self):
        return hash((self.x, self.y,self.z))
    
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
        self.center, self.radius = bound_sphere(v1,v2,v3)
        if not( np.isclose((self.center - v2).norm(),self.radius)):
            raise Exception('Error in bounding sphere calculation')
    
class Mesh:
    'Mesh class'
    def __init__(self):
        self.triangles = dict()
        self.size = 0
        self.KDTree = None
    
    
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
    
    def make_tree(self, depth = 8):
        'make KDTree structure'
        self.KDTree = KDTree(self.triangles.values(), depth = depth)
    
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
            c = closest_pt_on_triangle(tri,a)
            if (c-a).norm() < (closest_pt-a).norm():
                closest_pt = c
        return closest_pt
    
    def closest_pt_on_mesh(self, a):
        '''
        Find the closest point on mesh using either a KDTree search or a bf search.
        If the KDTree is constructed, then the function uses KDTree search, else it 
        calls bf_closet_pt_on_mesh
        Parameters
        ----------
        a : Vec3D
            Point in space

        Returns
        -------
        The closest point on mesh

        '''
        if self.KDTree != None:
            p, r = self.KDTree.find(a)
            return p
        else:
            return self.bf_closet_pt_on_mesh(a)
        



def closest_pt_on_triangle(tri, point):
    '''
    Finds the closest point on a given triangle from a given point

    Parameters
    ----------
    triangle : Triangle
    point : Vec3D
    Returns
    -------
    result: Vec3D
        The closest pt on the triangle
    '''
    p,r,q = tri.vertices
    #Configure lstsq
    A = np.hstack([(q-p).matrix, (r-p).matrix])
    b = (point-p).matrix
    res = np.linalg.lstsq(A,b,rcond=None)
    #unpack lambda and mu
    lam = res[0][0][0]
    miu = res[0][1][0]
    c = p+lam*(q-p)+miu*(r-p)
    #Check condition
    if lam >= 0 and miu >= 0 and (lam+miu) <= 1:
        pass
    else:
        #closest pt on edge
        if lam < 0:
            c = ProjectOnSegment(c,r,p)
        elif miu < 0:
            c = ProjectOnSegment(c,p,q)
        elif lam+miu > 1:
            c = ProjectOnSegment(c,q,r)
        
    return c

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
        
    
    
def bound_sphere(v1,v2,v3):
    '''
    Finds the smallest bounding sphere given three points in space

    Parameters
    ----------
    v1 : Vec3D
        first point
    v2 : Vec3D
        second point.
    v3 : Vec3D
        third point.

    Returns
    -------
    center : TYPE
        DESCRIPTION.
    r : TYPE
        DESCRIPTION.

    '''
    
    long_edge = sorted([(v1,v2),(v2,v3), (v1,v3)], key = lambda x:(x[0]-x[1]).norm())[0]
    a,b = long_edge
    c = (set([v1,v2,v3]) - set([a,b])).pop()
    #DEBUG
    if c == a or c == b:
        raise Exception()
    q = (a+b)/2
    #Try q with inequality check
    if (b-q).dot((b-q)) == (a-q).dot(a-q) and \
        (c-q).dot(c-q)<= (a-q).dot(a-q) and\
        ((b-a).cross(c-a)).dot(q-a) == 0:
            q = (a+b)/2
    else:
        f = (a+b)/2
        u = a-f
        v = c-f
        d = (u.cross(v)).cross(u)
        gama = (v.dot(v) - u.dot(u))/( (2*d).dot(v-u) )
        if gama <= 0:
            lam =0
        else:
            lam = gama
        q = f+lam*d
    r = (q-a).norm()
    return (q, r)
    
    

        

    

class KDTree:
    '3DTree'
    def __init__(self, List, depth = 8):
        x_sorted = sorted(List, key = lambda i : i.center.x)
        self.Lower = KDTreeNode(x_sorted[ : len(x_sorted)//2], 1 ,depth)
        self.Higher = KDTreeNode(x_sorted[len(x_sorted)//2 : ], 1, depth)
        self.midvalue = x_sorted[len(x_sorted)//2].center.x
        self.Lower_HighBound = self.midvalue + x_sorted[len(x_sorted)//2 - 1].radius
        self.Higher_LowerBound = self.midvalue - x_sorted[len(x_sorted)//2].radius
    
    def find(self,point):
        if point.x < self.midvalue:
            result  = self.Lower.find(point)
            (cur_opt_c,cur_opt_distance) = result
            if point.x + cur_opt_distance > self.Higher_LowerBound:
                return min(result, self.Higher.find(point),key = lambda x:x[1])
        else:
            result  = self.Higher.find(point)
            (cur_opt_c,cur_opt_distance) = result
            if point.x - cur_opt_distance < self.Lower_HighBound:
                return min(result, self.Lower.find(point),key = lambda x:x[1])
        
        return result
    
        
    
    
    
class KDTreeNode:
    'KDTree Node Data Stracture'
    def __init__(self,List,level,depth):
        #Check if is leaf node
        self.isLeaf = (level == depth-1) or len(List) == 1
        #determine which axis to split
        self.level = level
        self.split_by = 'xyz'[level % 3]
        if self.split_by == 'x':
            sort_key = lambda i:i.center.x
        if self.split_by == 'y':
            sort_key = lambda i:i.center.y
        if self.split_by == 'z':
            sort_key = lambda i:i.center.z
        sL = sorted(List, key = sort_key)
        #If is leaf node then store rest of list
        if self.isLeaf:
            self.container = List
        else:
            #Construct sub trees
            self.Lower = KDTreeNode(sL[: len(sL)//2], level+1, depth)
            self.Higher = KDTreeNode(sL[len(sL)//2 : ], level+1, depth)
            #Detemine Mid value
            if self.split_by == 'x':
                self.midvalue = sL[len(sL)//2].center.x
            if self.split_by == 'y':
                self.midvalue = sL[len(sL)//2].center.y
            if self.split_by == 'z':
                self.midvalue = sL[len(sL)//2].center.z
            #Determine subtree bounding box
            self.Lower_HighBound = self.midvalue + sL[len(sL)//2 - 1].radius
            self.Higher_LowerBound = self.midvalue - sL[len(sL)//2].radius
    
    
    def find(self,point):
        'Find the closest point on triangle'
        #If leaf node, calc potential closest point
        if self.isLeaf:
            all_closest_pts = sorted([ closest_pt_on_triangle(tri, point) for tri in self.container],key = lambda x:(point-x).norm())
            cur_opt_c = all_closest_pts[0]
            cur_opt_distance = (point-cur_opt_c).norm()
            return (cur_opt_c,cur_opt_distance)
            
            
        #Recursive finding leaf node and check finding returned by child
        if self.split_by == 'x':
            if point.x < self.midvalue:
                result  = self.Lower.find(point)
                (cur_opt_c,cur_opt_distance) = result
                if point.x + cur_opt_distance > self.Higher_LowerBound:
                    return min(result, self.Higher.find(point),key = lambda x:x[1])
            else:
                result  = self.Higher.find(point)
                (cur_opt_c,cur_opt_distance) = result
                if point.x - cur_opt_distance < self.Lower_HighBound:
                    return min(result, self.Lower.find(point),key = lambda x:x[1])
        if self.split_by == 'y':
            if point.y < self.midvalue:
                result  = self.Lower.find(point)
                (cur_opt_c,cur_opt_distance) = result
                if point.y + cur_opt_distance > self.Higher_LowerBound:
                    return min(result, self.Higher.find(point),key = lambda x:x[1])
            else:
                result  = self.Higher.find(point)
                (cur_opt_c,cur_opt_distance) = result
                if point.y - cur_opt_distance < self.Lower_HighBound:
                    return min(result, self.Lower.find(point),key = lambda x:x[1])
        if self.split_by == 'z':
            if point.z < self.midvalue:
                result  = self.Lower.find(point)
                (cur_opt_c,cur_opt_distance) = result
                if point.z + cur_opt_distance > self.Higher_LowerBound:
                    return min(result, self.Higher.find(point),key = lambda x:x[1])
            else:
                result  = self.Higher.find(point)
                (cur_opt_c,cur_opt_distance) = result
                if point.z - cur_opt_distance < self.Lower_HighBound:
                    return min(result, self.Lower.find(point),key = lambda x:x[1])
    
        return result
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    