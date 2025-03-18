"""
Mol2Mesh
Main function

Written by Ira Novianti
"""

import itertools
import numpy as np
import webcolors
import plotly.graph_objects as go
from sympy import Plane, Point3D, Line3D

from .file_parser import *
from .style import *
from .surfaces import *

class Mol2Mesh:
	def __init__(self, file_path:str, style='BallStick', multicov=False, res_a:int=25, res_b:int=15, name:str=None):
		"""
		A Mol2Mesh object contains triangular meshes of spheres (atoms) and cylinders (bonds)

		Parameters
		file_path : str
			path to chemical table file (sdf, mol, or gjf)
		style: str or class
			style configuration defining radii and colors of atoms and bonds
		res_a: int
			resolution of an atom sphere (number of outer points in each vertical and horizontal circles defining a sphere)
		res_b: int
			resolution of a bond cylinder (number of points of each top and bottom bases of a cylinder)
		name : str
			the name of the molecule which will be stored in the resulting 3D model stl file
		"""

		self.name = name
		if name is None:
			self.name = file_path.split('/')[-1].split('.')[0]
		
		#1. read chemical table file
		parser = {'mol': open_mol, 'sdf': open_mol, 'gjf': open_gjf}

		file_type = file_path.split('.')[-1]
		assert file_type in parser, "File type is not supported yet"

		self.mol_d = parser[file_type](file_path)
		
		#2. define style
		style_d = {'BallStick': BallStick(), 'Sticks': Sticks(), 'SpaceFilling': SpaceFilling()}

		style_is_dict = style in style_d
		style_is_custom = hasattr(style, '__dict__')

		assert style_is_dict or style_is_custom,\
		"Style should be a configuration class or a string of one of the defined styles ['BallStick', 'Sticks', 'SpaceFilling']"

		self.style = style_d[style] if style_is_dict else style

		#3. calculate triangular surface of spheres of atoms and cylinders of bonds
		for atom in self.mol_d['atoms']:
			element = atom['element'] if atom['element'] in self.style.atom_radius else 'others'
			R = self.style.atom_radius[element]
			verts, fcs = sphere_tri(R=R, res=res_a, origin=atom['coor'])

			element = atom['element'] if atom['element'] in self.style.atom_color else 'others'
			color = self.style.atom_color[element]
			atom.update({'vertices': verts,
						 'faces': fcs,
						 'color': tuple(webcolors.name_to_rgb(color))})
		
		for bond in self.mol_d['bonds']:

			R = self.style.BOND_RADIUS

			if multicov and bond.get('bond_type', None) in ['double', 'triple']:

				if bond['bond_type']=='double':
					_R = 0.47 * R
					dist = 0.27 * R
				elif bond['bond_type']=='triple':
					_R = 0.3 * R
					dist =  R

				bond_list = [[tuple(b['coor_1']), tuple(b['coor_2'])] for b in self.mol_d["bonds"]]
				multi_bond = multiple_bond(bond['coor_1'], bond['coor_2'], bond_list, dist=dist, bond_type=bond['bond_type'])

				verts = []
				fcs = []

				for sub_bond in multi_bond:
					sub_verts, sub_fcs = cylinder_tri(*sub_bond, R=_R, res=res_b)

					init_idx = len(verts)
					verts.extend(sub_verts)
					fcs.extend(sub_fcs+init_idx)
				
				verts = np.array(verts)
				fcs = np.array(fcs)

			else:
				verts, fcs = cylinder_tri(bond['coor_1'], bond['coor_2'], R=R, res=res_b)
			
			color = self.style.BOND_COLOR
			bond.update({'vertices': verts,
						 'faces': fcs,
						 'color': tuple(webcolors.name_to_rgb(color))})
		
		self.signature = "created by Mol2Mesh (c) 2024 github.com/iranovianti"
	
	def get_mesh(self, as_type='verts_and_faces'):
		"""
		Get a combined array of meshes
		If as_type: 'verts_and_faces', returns a combined arrays of vertices and faces
		(vertices are in (x, y, z) coordinates, faces are vertex indexes defining triangles)

		If as_type: 'triangles', returns a combined arrays of triangles
		(each triangle consists of 3 points (p1, p2, p3), p is xyz coodinates)
		triangles = vertices[faces]
		"""
		assert as_type in ['verts_and_faces', 'triangles']

		vertices = []
		faces = []

		for surf in itertools.chain(self.mol_d['atoms'], self.mol_d['bonds']):
			init_idx = len(vertices)
			vertices.extend(surf['vertices'])
			faces.extend(surf['faces'] + init_idx)
		
		vertices = np.array(vertices)
		faces = np.array(faces)

		if as_type=='verts_and_faces':
			return vertices, faces
		
		elif as_type=='triangles':
			return vertices[faces]
	
	def show(self):
		"""
		3D visualization of the molecule with plotly
		"""
		fig = go.Figure()
		for surf in itertools.chain(self.mol_d['atoms'], self.mol_d['bonds']):
			x,y,z = surf['vertices'].T
			i,j,k = surf['faces'].T
			fig.add_trace(go.Mesh3d(x=x, y=y, z=z,
									i=i, j=j, k=k,
									color=f"rgb{surf['color']}"))
		fig.update_scenes(aspectmode='data', xaxis_visible=False, yaxis_visible=False, zaxis_visible=False)
		fig.show()
	
	def show_mesh(self, show_edges=True, show_vertices=False):
		"""
		Show molecule meshes (vertices and edges) as wireframes by plotly
		based on https://community.plotly.com/t/show-edges-of-the-mesh-in-a-mesh3d-plot/33614/3
		"""
		vertices, faces = self.get_mesh()

		triangles = vertices[faces]

		fig = go.Figure()

		if show_edges:
			Xe = []
			Ye = []
			Ze = []
			for T in triangles:
				Xe.extend([T[k%3][0] for k in range(4)]+[ None])
				Ye.extend([T[k%3][1] for k in range(4)]+[ None])
				Ze.extend([T[k%3][2] for k in range(4)]+[ None])
			
			fig.add_trace(go.Scatter3d(x=Xe,y=Ye,z=Ze,
									   mode='lines',
									   line=dict(color='gray')))
		
		if show_vertices:
			fig.add_trace(go.Scatter3d(x=vertices[:,0],
									   y=vertices[:,1],
									   z=vertices[:,2],
									   mode='markers',
									   marker=dict(size=2, color='black')))
		
		fig.update_scenes(aspectmode='data', xaxis_visible=False, yaxis_visible=False, zaxis_visible=False)
		fig.show()

	
	def save_mesh(self, filename):
		"""
		Save molecule mesh as an .stl file.
		The resulting file can be used for 3D printing.
		stl files do not store color. If you want a 3D model with colors, try save_glb()
		"""
		if filename.split('.')[-1]!='stl':
			filename = filename+'.stl'
		
		f = open(filename, 'w')
		f.write(f'solid {self.name} - {self.signature}\n')

		triangles = self.get_mesh(as_type='triangles')

		for p1, p2, p3 in triangles:
			if np.isnan(np.array([p1,p2,p3])).any():
				continue
			n = find_normal(p1, p2, p3)
			f.write(f'facet normal {n[0]:.7f} {n[1]:.7f} {n[2]:.7f}\n')
			f.write('outer loop\n')
			f.write(f'vertex {p1[0]:.7f} {p1[1]:.7f} {p1[2]:.7f}\n')
			f.write(f'vertex {p2[0]:.7f} {p2[1]:.7f} {p2[2]:.7f}\n')
			f.write(f'vertex {p3[0]:.7f} {p3[1]:.7f} {p3[2]:.7f}\n')
			f.write('endloop\n')
			f.write('endfacet\n')
		
		f.write(f'endsolid {self.name} - {self.signature}\n')
		f.close()
	
	def save_glb(self, filename):
		"""
		Save the molecule as a 3D model with colors (.glb).
		This function relies on trimesh library
		"""
		if filename.split('.')[-1]!='glb':
			filename = filename+'.glb'
		
		import trimesh

		meshes = []

		for surf in itertools.chain(self.mol_d['atoms'], self.mol_d['bonds']):
			mesh = trimesh.Trimesh(vertices=surf['vertices'], faces=surf['faces'])
			mesh.visual.face_colors = np.array((*surf['color'], 255))
			meshes.append(mesh)
		
		combined = trimesh.util.concatenate(meshes)
		combined.export(filename)

#---Helper functions for creating multiple covalent bonds---

def get_reference_atom(A, B, bond_list):
	"""
	get a reference atom to form a plane for double bond (sp2 hybridization)
	by finding a neighboring atom with the most covalent bonds
	"""
	#sort atoms based on the number of covalent bonds
	bonds_flat = list(itertools.chain.from_iterable(bond_list))
	atms_sorted = sorted(set(bonds_flat), key = lambda atm: bonds_flat.count(atm))[::-1]

	A = tuple(A)
	B = tuple(B)

	#atoms bonded to atom A and B
	A_bonds = [a for b in bond_list if A in b and B not in b for a in b if a!=A]
	B_bonds = [a for b in bond_list if B in b and A not in b for a in b if a!=B]

	#get reference atom: an atom bonded to either A or B with the most covalent bonds
	C = None

	for atm in atms_sorted:
		if atm in A_bonds + B_bonds:
			C = atm
			break
	return C

def multiple_bond(A, B, bond_list, dist=0.1, bond_type='double'):
	"""
	turn bond coordinates into multiple coordinates of sub bonds
	"""
	if bond_type=='double':
		"""
		Create a perpendicular plane (`plane2`) by using the original bond as the normal,
		and make the plane intersection be the direction in which the coordinates are moved
		to make sub bonds.

		-contributed by https://github.com/mziele1
		"""
		#form sp2 plane
		C = get_reference_atom(A, B, bond_list)
	
		pA = Point3D(A)
		pB = Point3D(B)
		pC = Point3D(C)

		plane = Plane(pA, pB, pC) #the plane of target bond

		#get a plane perpendicular to sp2 plane
		normal = Line3D(pA, pB)
		plane2 = Plane(pA, normal_vector=normal.direction)
	
		#get direction to multiply target bond
		direction = plane.intersection(plane2)[0].direction
	
	elif bond_type=='triple':
		#get random point perpendicular to the bond
		pts = dir_disk(A, B, R=1., res=10)
		C = [pts[0][i] + A[i] for i in [0, 1, 2]]
		#get direction to multiply target bond
		direction = Line3D(A, C).direction
		
	dir_vec = np.array([float(direction.x),
						float(direction.y),
						float(direction.z)])
	
	A1 = A + dir_vec * dist
	B1 = B + dir_vec * dist

	A2 = A - dir_vec * dist
	B2 = B - dir_vec * dist

	if bond_type == 'double':
		return (A1, B1), (A2, B2)
    
	elif bond_type == 'triple':
		return (A, B), (A1, B1), (A2, B2)