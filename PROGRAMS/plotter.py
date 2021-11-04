# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 04:59:41 2021

@author: MikeH
"""

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt


fig = plt.figure()
ax = plt.axes(projection='3d')

def plot_vec_list(vec_list, color='b'):
    ax.scatter([p.x for p in vec_list], [p.y for p in vec_list], [p.z for p in vec_list], c = color)

def plot_arrow(start, vec):
    A = np.hstack([start.matrix.T, vec.matrix.T])[0]
    print(A)
    ax.quiver(A[0],A[1],A[2],A[3],A[4],A[5])