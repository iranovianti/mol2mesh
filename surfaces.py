"""
Mol2Mesh
Functions to make triangular surfaces

Written by Ira Novianti
"""

import numpy as np
from scipy.spatial import Delaunay
from scipy.linalg import norm
import math

def sphere_tri(R=1., res=10, origin=[0, 0, 0]):
	phi = np.linspace(0, 2*np.pi, res)
	theta = np.linspace(0, np.pi, res)

	[X, Y] = np.meshgrid(phi, theta)
	#perform triangulation
	triangles=Delaunay(np.array([X.flatten(), Y.flatten()]).T).simplices

	#sphere vertices
	x = R * np.outer(np.cos(phi), np.sin(theta))
	y = R * np.outer(np.sin(phi), np.sin(theta))
	z = R * np.outer(np.ones_like(phi), np.cos(theta))

	#translate sphere
	x = (x+origin[0]).flatten()
	y = (y+origin[1]).flatten()
	z = (z+origin[2]).flatten()

	points = np.array([x,y,z]).T

	return points, triangles


def cylinder_tri(p0, p1, R=1., res=10):
	disk = dir_disk(p0, p1, R, res)
	disk0 = np.array([p0[i] + disk[i] for i in [0, 1, 2]]) #disk translated to the top base
	disk1 = np.array([p1[i] + disk[i] for i in [0, 1, 2]]) #disk translated to the bottom base

	#combine all vertices
	points = np.concatenate([[p0], disk0.T, [p1], disk1.T])

	tri0 = [[0, i, i+1] for i in range(1,res)] #triangles for top base
	tri1 = np.array(tri0) + res + 1 #triangles for bottom base
	tri2 = [] #triangles for the side
	for i in range(1, res):
		j = i + res + 1
		tri2.extend([[i, j, i+1], [i+1, j, j+1]])

	#combine all triangles
	triangles = np.concatenate([tri0, tri1, tri2])

	return points, triangles

def dir_disk(p0, p1, R, res):
	p0 = np.array(p0)
	p1 = np.array(p1)

	v = p1 - p0
	v = v / norm(v)

	not_v = np.array([1, 0, 0])
	if (np.absolute(v) == not_v).all():
		not_v = np.array([0, 1, 0])

	n1 = np.cross(v, not_v)
	n1 /= norm(n1)

	n2 = np.cross(v, n1)

	theta = np.linspace(0, 2*np.pi, res)

	disk = [R * np.sin(theta) * n1[i] + R * np.cos(theta) * n2[i] for i in [0, 1, 2]]
	return disk

def find_normal(p1, p2, p3):
	v1 = p2 - p1
	v2 = p3 - p1
	v3 = np.cross(v1, v2)
	n = v3 / math.sqrt(np.sum(v3*v3))
	return n
