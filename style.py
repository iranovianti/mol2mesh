"""
Mol2Mesh
Base style configuration

Written by Ira Novianti
"""

class style_config:
	ATOM_SCALE = 0.25
	ATOM_RADII = None
	BOND_RADIUS = 0.2
	Rvdw = {'C': 1.7, 'H': 1.2, 'O': 1.55, 'N': 1.6, 'P': 1.95, 'S': 1.8,
			'Ca': 2.4, 'Mg': 2.2, 'Li': 2.2, 'Na': 2.4, 'K': 2.8,
			'F': 1.5, 'Cl': 1.8, 'Br': 1.9, 'I': 2.1,
			'He': 1.4, 'Ne': 1.54, 'Ar': 1.88, 'Kr': 2.02, 'Xe': 2.16,
			'Zn': 2.1, 'Si': 2.1, 'B': 1.8, 'Al': 2.1, 'Ni': 2.0, 'Cu': 2.0, 'Se': 1.9, 'As': 2.05,
			'Ag': 2.1, 'Au': 2.1, 'Pd': 2.05, 'Pt': 2.05, 'Pb': 2.3, 'Hg': 2.05, 'others': 2.0}
	
	colors = {'N': 'blue',
			'C': 'gray',
			'O': 'red',
			'H': 'white',
			'Cl': 'yellow',
			'others': 'black'}
	
	BOND_COLOR = 'gray'
	ATOM_COLORS = None
	
	def __init__(self):
		if self.ATOM_RADII:
			self.atom_radius = {k: self.ATOM_RADII for k in self.Rvdw.keys()}
		else:
			self.atom_radius = {k: self.ATOM_SCALE * v for k,v in self.Rvdw.items()}
		
		if self.ATOM_COLORS:
			self.atom_color = {k: self.ATOM_COLORS for k in self.colors.keys()}
		else:
			self.atom_color = self.colors

class BallStick(style_config):
	pass

class Sticks(style_config):
	RADIUS = 0.2   
	BOND_RADIUS = RADIUS
	ATOM_RADII = RADIUS

	COLOR = 'gray'
	BOND_COLOR = COLOR
	ATOM_COLORS = COLOR


class SpaceFilling(style_config):
	ATOM_SCALE = 1.