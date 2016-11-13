#!/usr/bin/env python

import sys
import numpy as np
import argparse as arg
from fractions import gcd
from matplotlib import path
from matplotlib import gridspec
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle, PathPatch


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='''
    Generates a Grid of points according to two basis vectors and rolls it
    according to a Chiral Vector.''')

    # Optional arguments
    parser.add_argument('-a1', '--a1', default=15.6, type=float, help='''
    Length of Basis Vector 1.''')

    parser.add_argument('-a2', '--a2', default=9.5, type=float, help='''
    Length of Basis Vector 2.''')

    parser.add_argument('-g', '--gamma', default=90.0, type=float, help='''
    Angle between the two Basis Vectors.''')

    parser.add_argument('-n', '--n', default=2, type=int, help='''
    Chiral Index of Basis Vector 1.''')

    parser.add_argument('-m', '--m', default=3, type=int, help='''
    Chiral Index of Basis Vector 2.''')

    parser.add_argument('-r', '--rep', default=1, type=int, help='''
    Number of repetitions of the Tubular unit cell.''')

    parser.add_argument('-o', '--output', default='tube', help='''
    Root of the name of the output file.''')

    args = parser.parse_args()

    return args


def v1v2_angle(v1, v2):
    '''Returns the angle between two vectors.'''
    # Remember that the angle between a plane and a vector equals
    # 90 - alpha, where alpha is the angle between the vector and
    # the normal to the plane. To obtain such an angle, you could
    # do angle = 90 - v1v2_angle(v1, np.cross(x, y)), where x, y
    # are the two vectors that define the plane.

    dotprod = np.dot(v1, v2)
    try:
        theta = np.degrees(np.arccos(dotprod / (np.linalg.norm(v1) * np.linalg.norm(v2))))
    except:
        theta = 0.0

    return theta


def rot(axis, theta):
    '''Returns the rotation matrix for the anticlockwise rotation about
    axis by theta according to Rodrigues' formula.'''

    axis = axis / np.linalg.norm(axis)
    theta = -1 * np.radians(theta)
    I = np.eye(3)

    # Define axis' cross-product matrix
    K = np.cross(I, axis)

    R = I + np.sin(theta) * K + (1 - np.cos(theta)) * np.linalg.matrix_power(K, 2)

    return R


if __name__ == '__main__':

    args = options()
    origin = np.zeros(3)
    cartesian = np.eye(2)
    ux = cartesian[0]
    uy = cartesian[1]
    reps = args.rep

    #
    # Parallelogram data in Angstrom
    #
    side1 = args.a1
    side2 = args.a2
    gamma = np.radians(args.gamma)
    n = args.n
    m = args.m
    
    a1 = ux * side1
    a2 = np.array([side2 * np.cos(gamma), side2 * np.sin(gamma)])
    
    #
    # Tubular cell
    # C: Chiral Vector
    # T: Translation Vector
    #
    dR = gcd( (2 * m + n), (m + 2 * n))
    C = n * a1 + m * a2
    t1 = (2 * m + n) / dR
    t2 = (2 * n + m) / dR
    T = t1 * a1 - t2 * a2

    r = np.linalg.norm(C) / (2 * np.pi)
    
    #
    # Define Edges of the Tubular Cell
    #
    O = origin[:-1]
    A = C
    B = C + T
    D = T
    rectangle = path.Path([O, A, B, D])
    
    #
    # Generate points inside the Tubular cell
    #
    points = []
    for i in range(n * t1):
        for j in range(m * t2):
    
            x = i * side1 + j * side2 * np.cos(gamma)
            y = j * side2 * np.sin(gamma)
            point = np.array([x, y])
    
            if rectangle.contains_points(np.atleast_2d(point)):
                points.append(point)
    
    
            xneg = i * side1 - j * side2 * np.cos(gamma)
            yneg = -j * side2 * np.sin(gamma)
            point = np.array([xneg, yneg])
    
            if rectangle.contains_points(np.atleast_2d(point)):
                points.append(point)
    
    points = np.array(points)

    xs = points[:,0]
    ys = points[:,1]
    
    # fig = plt.figure(figsize=(16, 12))
    # gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1]) 
    # ax0 = plt.subplot(gs[0])
    # ax0.scatter(xs,ys, color="k", s=1)
    # ax0.plot([D[0], B[0]], [D[1], B[1]], color="b", lw=2, label="Chiral Vector")
    # ax0.plot([O[0], A[0]], [O[1], A[1]], color="b", lw=2)
    # ax0.plot([D[0], O[0]], [D[1], O[1]], color="g", lw=2, label="Translation Vector")
    # ax0.plot([B[0], A[0]], [B[1], A[1]], color="g", lw=2)
    # ax0.legend().draw_frame(False)
    # plt.show()

    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1]) 
    
    #
    # Rotate everything to superimpose the Chiral Vector with the x axis
    #
    angle = v1v2_angle(C, ux)
    uz = np.array([0, 0, 1])
    M = rot(uz, -angle)[:-1,:-1]
    
    #
    # Transformation to have the lowest left edge of the tubular cell correspond
    # to the origin
    #
    D = np.dot(D, M)
    O = np.dot(O, M) - D
    A = np.dot(A, M) - D
    B = np.dot(B, M) - D
    points = np.dot(points, M) - D
    D -= D
    xs = points[:,0]
    ys = points[:,1]
    
    ax0 = plt.subplot(gs[0])
    ax0.scatter(xs,ys, color="k", s=1)
    ax0.plot([D[0], B[0]], [D[1], B[1]], color="b", lw=2, label="Chiral Vector")
    ax0.plot([O[0], A[0]], [O[1], A[1]], color="b", lw=2)
    ax0.plot([D[0], O[0]], [D[1], O[1]], color="g", lw=2, label="Translation Vector")
    ax0.plot([B[0], A[0]], [B[1], A[1]], color="g", lw=2)
    ax0.plot([D[0], D[0]], [D[1], O[1]], color="k", lw=2, label="Tube Axis")
    ax0.legend().draw_frame(False)

    #
    # Height of the Tube for repetitions
    #
    h = O[1] - D[1]
    
    #
    # Roll Tubular cell
    #
    rolled = []
    for point in points:
    
        fracx = point[0] / B[0]
        phi = 2 * np.pi * fracx
        rolledx = r * np.cos(phi)
        rolledy = r * np.sin(phi)
        rolledz = point[1]
    
        rolledpoint = np.array([rolledx, rolledy, rolledz])
        rolled.append(rolledpoint)
    
    rolled = np.array(rolled)

    for rep in range(reps):

        z = rolled[:,2] + rep * h
        repeated = np.copy(rolled)
        repeated[:,2] = z
        rolled = np.r_[ rolled, repeated ]
   
    x = rolled[:,0]
    y = rolled[:,1]
    z = rolled[:,2]

    ax1 = plt.subplot(gs[1], projection= '3d', aspect='equal')
    ax1.scatter(x,y,z, marker=".", color="k")
    ax1.set_xlim(x.min() - 1, x.max() + 1)
    ax1.set_ylim(y.min() - 1, y.max() + 1)
    ax1.set_zlim(z.min() - 1, z.max() + 1)
    ax1.plot([0, 0], [0, 0], [0, z.max()], color="k", lw=2)
    circle = Circle((0, 0), r, color="b", lw=2, fill=False)
    ax1.add_patch(circle)
    art3d.pathpatch_2d_to_3d(circle, z=-3, zdir="z")

    atoms = np.ones(len(rolled))
    rolled = np.c_[atoms, rolled]

    #
    # Sort along z for easy cutting
    #
    rolled = rolled[rolled[:, -1].argsort()]

    with open(args.output + ".coords.xyz", "w") as f:
        f.write("%d\n\n" % len(rolled))
        np.savetxt(f, rolled, fmt="%5d %14.8f %14.8f %14.8f")

    # plt.show()
