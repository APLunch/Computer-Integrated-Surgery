# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 12:59:30 2021

@author: Yunxin Chen
"""

import numpy as np
import cismath as cis 
from registration import registration
import random
import math
from scipy.optimize import least_squares
np.set_printoptions(precision=3)

# read the calbody txt file
def read_body(filename):
    
    with open("../Input Data/{}".format(filename),'r') as f:
        lines = f.readlines()
        
        #read first line to get ND, NA and NC
        first_line = lines[0].split(" ")
        N_markers = int(first_line[0])
        
        #save all the coordinates of d, a, c to three defined lists of Vec3D object
        points_LED = []
        
        for line in lines[1:N_markers+1]:
            words = line.split()
            words = [word.strip(' ') for word in words]
            #print(float(words[0]), float(words[1]), float(words[2]))
            points_LED.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
        
        words = lines[N_markers+1].split()
        words = [word.strip(' ') for word in words]
        #print(float(words[0]), float(words[1]), float(words[2]))
        p_tip = cis.Vec3D(float(words[0]), float(words[1]), float(words[2]))
    
    return N_markers, points_LED, p_tip


# read the calreadings txt file
def read_sample_readings(filename, N_A, N_B):
    
    with open("../Input Data/{}".format(filename),'r') as f:
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
        
                #print(float(words[0]), float(words[1]), float(words[2]))
        
            for line in lines[1+N_A+i*(N_S) : N_A+1+N_B+i*(N_S)]:
                words = line.split()
                words = [word.strip(',') for word in words]
                B_list.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
                
                #print(float(words[0]), float(words[1]), float(words[2]))
                
            for line in lines[N_A+1+N_B+i*(N_S) : N_A+1+N_B+N_D+i*(N_S)]:
                words = line.split()
                words = [word.strip(',') for word in words]
                D_list.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
        
                #print(float(words[0]), float(words[1]), float(words[2]))
                
    return A_list, B_list, D_list, N_D, N_samples

def load_mesh_from_file(filename):
    Mesh = cis.Mesh()
    with open("../Input Data/{}".format(filename),'r') as f:
        lines = f.readlines()
        #Get number of vertices
        N_Vertex = int(lines[0].strip())
        #Get vertices
        vertices = []
        for i in range(1,N_Vertex+1):
            line = lines[i].strip()
            coordinate = [float(word.strip()) for word in line.split()]
            vertices.append(cis.Vec3D(coordinate[0],coordinate[1],coordinate[2]))
        Mesh.vertices = vertices
        #Get number of triangles
        N_triangles = int(lines[N_Vertex+1].strip())
        #Get triangles
        for i in range(N_Vertex+2,len(lines)):
            line = lines[i].strip()
            #get index of vertices
            i1,i2,i3 = [int(word.strip()) for word in lines[i].split()[:3] ]
            v1,v2,v3 = [vertices[index] for index in [i1,i2,i3]]
            tri = cis.Triangle(v1,v2,v3, i1,i2,i3)
            Mesh.add(tri)
    return Mesh



# This function read the C value from the given output txt file
def get_debug_output(file_name):
    #read txt file line by line
    lines = []
    with open("../Input Data/{}".format(file_name),'r') as f:
        lines = f.readlines()
    f.close()
    
    #read first line to get N_C and N_frames
    first_line = lines[0].split(" ")
    N_samps = int(first_line[0])
    
    d = []
    c = []
    error = []
    for line in lines[1:N_samps+1]:
        words = line.split()
        words = [word.strip(',') for word in words]
        #print(float(words[0]), float(words[1]), float(words[2]))
        d.append(cis.Vec3D(float(words[0]), float(words[1]), float(words[2])))
        c.append(cis.Vec3D(float(words[3]), float(words[4]), float(words[5])))
        error.append(float(words[6]))

    return d, c, error


def compare_output(d_list, c_list, e, d_output, c_output, e_output, N_samples):
    e = []
    for row in range(N_samples):
        e.append(abs((d_list[row] - c_list[row]).norm()))
    
    sum_error_d = cis.Vec3D(0, 0, 0)
    sum_error_c = cis.Vec3D(0, 0, 0)
    sum_error_e = 0
    for row in range(N_samples):
        sum_error_d += cis.Vec3D(abs(d_list[row].x - d_output[row].x), abs(d_list[row].y - d_output[row].y), abs(d_list[row].z - d_output[row].z))
        sum_error_c += cis.Vec3D(abs(c_list[row].x - c_output[row].x), abs(c_list[row].y - c_output[row].y), abs(c_list[row].z - c_output[row].z))
        sum_error_e += abs(e[row] - e_output[row])
        
    average_error_d = cis.Vec3D(sum_error_d.x / N_samples, sum_error_d.y / N_samples, sum_error_d.z / N_samples)
    average_error_c = cis.Vec3D(sum_error_c.x / N_samples, sum_error_c.y / N_samples, sum_error_c.z / N_samples)
    average_error_e = sum_error_e / N_samples

    return average_error_d, average_error_c, average_error_e

         
        
        
def ICP(d_list, Mesh):
    # Step 0: Initialization
    F = []
    # initialize F0
    F.append(cis.Frame(cis.Rot3D(np.eye(3)),cis.Vec3D(0,0,0)))
    threshold = []
    # initialize threshold0
    threshold.append(5.0)
    sigma = []
    #sigma.append(0.0)
    error_max = []
    #error_max.append(0.0)
    error_mean = []
    #error_mean.append(0.0)
    
    check = 0
    n = 0
    d_list_sample = random.sample(d_list,len(d_list)//4)
    
    
    #Input: M, d, F0, n0 check
    iteration_with_sample = 0
    while check < 8:
        # step 1: matching
        # create two empty lists A and B, which is different from Body A and Body B
        A_group = []
        B_group = []
        c_list = []
        s_list = []
        e_list = []
        for d in d_list_sample:
            s = F[n]*d
            s_list.append(s)
            c,tri = Mesh.closest_pt_on_mesh(s)
            c_list.append(c)
            e = (s -c).norm()
            e_list.append(e)
            if e < threshold[n]:
                A_group.append(d)
                B_group.append(c)
        
        # step 2: transformation part
        F.append(registration(A_group, B_group))
        
        E_list = []
        for k in range(len(A_group)):
            E = B_group[k] - F[n] * A_group[k]
            E_list.append(E)
        
        
        #calculate sigma, error_max and mean error
        sum_dot = 0
        sum_root = 0
        em = 0 
        for k in range(len(E_list)):
            e_dot = E_list[k].dot(E_list[k])
                
            sum_dot += e_dot
            
            if em < math.sqrt(e_dot):
                em = math.sqrt(e_dot)
                
            sum_root += math.sqrt(e_dot)
            
        sigma.append(math.sqrt(sum_dot)/len(E_list))
        error_max.append(em)
        error_mean.append(sum_root/len(E_list))
        
        
        
        # Step 3: adjustment
        #adjust threshold
        eta = 3*error_mean[n]
        threshold.append(eta)
        #print(error_mean[n])
        
        # Step 4 :(iteration) Termination 
        if 0.98 <= error_mean[n]/error_mean[n-1] <= 1.005:
            check += 1
            if check == 7 and d_list_sample != d_list:
                d_list_sample = d_list
                iteration_with_sample = n
        else:
            check = 0
        
        #Sigma, error max and error mean termination condition
        n += 1
        if sigma[-1] < 0.0008 and error_mean[-1] < 0.0008 and error_max[-1] < 0.0008:
            if d_list_sample != d_list:
                d_list_sample = d_list
            else:
                break
        #Hard Change if not converging
        if n > 75 and d_list_sample != d_list:
            iteration_with_sample = n
            d_list_sample = d_list
        if n > 150:
            break
        
    print("Error_Max:",error_max[-1])
    print("Sigma:",sigma[-1])
    print("Error_Mean:",error_mean[-1])
    print("Iteration with sampled d: {}".format(iteration_with_sample))
    print("Total Iteration: {}".format(n))
    return (F[-1],s_list,c_list,e_list, error_max[-1],error_mean[-1])
        
def read_modes(fname):
    '''
    Read in Modes of the mesh deformation. Meant to be used in PA5

    Parameters
    ----------
    fname : str
        filename containning the modes infomation

    Returns
    -------
    list(list(cis.Vec3D))

    '''
    with open('../Input Data/{}'.format(fname)) as f:
        lines = f.readlines()
        modes = []
        #Read  first line
        N_vertex = (int)(lines[0].split()[1].split('=')[1].strip())
        N_modes = (int)(lines[0].split()[2].split('=')[1].strip())
        #Read modes
        for i in range(len(lines)):
            if 'Mode ' in lines[i]:
                mode = []
                for j in range(i+1, i+1+N_vertex):
                    x,y,z = [ (float)(word.strip(' ,')) for word in lines[j].split()]
                    mode.append(cis.Vec3D(x,y,z))
                modes.append(mode)
                i = i+N_vertex
                    
        
    return modes
        
        
def ICP_And_Calc_Deformation_Weights(d_list,Mesh,Modes):
    ICP_Error_Mean = 99999
    residual = 99999
    weights = np.array([76,-37,-9.98,159,-33,101])
    lowest_weights = None
    lowest_residual = 99999
    while ICP_Error_Mean > 0.001:
        F_reg, s_list,c_list,e_list,e_max, ICP_Error_Mean = ICP(d_list,Mesh)
        n = 0
        check = 0
        while check < 8: 
            c_tri_list = [ Mesh.closest_pt_on_mesh(s) for s in s_list ]
            c_list = [ x[0] for x in c_tri_list]
            tri_list = [ x[1] for x in c_tri_list]
            #lstsq calculate weights
            res = lstsq_calc_weights(s_list, c_list,tri_list, Modes)
            weights = res[0].reshape(-1,)
            new_residual = np.linalg.norm(res[1])
            #record lowest redisual
            if new_residual < lowest_residual:
                lowest_weights = weights
                lowest_residual = new_residual
            #update mesh
            Mesh.vertices = [vertex_with_modes(vi, Modes, weights) for vi in range(len(Mesh.vertices))]
            Mesh.update_triangles()
            Mesh.make_tree(depth = -1)
            #Termination condition
            n += 1
            print("========================")
            print("Weight:", weights)
            print("Current Best :",lowest_weights)
            print("Residual:",residual)
            print("========================")
            #Termination condition
            if new_residual/residual < 1.01 and new_residual/residual > 0.98:
                check += 1
            else:
                check = 0
            if new_residual < 1e-8:
                break
            residual = new_residual
            
            
    return weights
    


def lstsq_calc_weights(s_list, c_list, tri_list, Modes):
    q0_list = []
    #calculate barycetric coordinate for all point s
    bc_list = [ barycentric(c_list[i],tri_list[i]) for i in range(len(c_list)) ]
    #Construct q0 matrics
    for i,tri in enumerate(tri_list):
        bc = bc_list[i]
        v1,v2,v3 = [v for v in tri.v_index]#index of vertex of this triangle
        m1,m2,m3 = [ bc.x*Modes[0][v1], bc.y*Modes[0][v2], bc.z*Modes[0][v3]]
        q0_list.append( m1+m2+m3 )
    q0_mat = cis.vec_list_to_matrix(q0_list).T.reshape(-1,1)
    #Construct s matrix
    s_mat = cis.vec_list_to_matrix(s_list).T.reshape(-1,1)
    #Construct q matrix
    q_mat = None
    for i,tri in enumerate(tri_list):
        bc = bc_list[i]
        qi_moded=[]
        for j in range(1,len(Modes)):
            m1,m2,m3 = [Modes[j][v_index] for v_index in tri.v_index]
            q = bc.x*m1+bc.y*m2+bc.z*m3
            qi_moded.append(q)
        qi_mat = cis.vec_list_to_matrix(qi_moded)
        if q_mat is None:
            q_mat = qi_mat
        else:
            q_mat = np.vstack([q_mat,qi_mat])
    #print(s_mat-q0_mat)
    x = np.linalg.lstsq(q_mat,s_mat-q0_mat , rcond=-1 )
    return x

def barycentric(c,tri):
    m1,m2,m3 = [ v for v in tri.vertices ]
    m_matrix = cis.vec_list_to_matrix([m1,m2,m3])
    barycentric_coord = cis.Vec3D(np.matmul(np.linalg.inv(m_matrix), c.matrix))
    return barycentric_coord
    


def vertex_with_modes(v_index, modes, weights):
    '''
    Apply deformation modes to a single vertex
    
    Parameters
    ----------
    v_index: int
        Index of a vertex
    modes : list(list(cis.Vec3D))
        modes of deformation
    weights : list(float)
        weights of modes

    Returns
    -------
    v: cis.Vec3D
        The vertex vector after mmodes 
    '''
    v = modes[0][v_index]
    for i in range(1,1+len(weights)):
        v = v + (weights[i-1]*modes[i][v_index])
    return v
            
        
            
        
        
        