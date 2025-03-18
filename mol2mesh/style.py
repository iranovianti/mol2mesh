"""
Mol2Mesh
Base style configuration

Written by Ira Novianti
"""

class style_config:
	ATOM_SCALE = 0.25
	BOND_RADIUS = 0.2
	Rvdw = {'C': 1.7, 'H': 1.2, 'O': 1.55, 'N': 1.6, 'P': 1.95, 'S': 1.8,
			'Ca': 2.4, 'Mg': 2.2, 'Li': 2.2, 'Na': 2.4, 'K': 2.8,
			'F': 1.5, 'Cl': 1.8, 'Br': 1.9, 'I': 2.1,
			'He': 1.4, 'Ne': 1.54, 'Ar': 1.88, 'Kr': 2.02, 'Xe': 2.16,
			'Zn': 2.1, 'Si': 2.1, 'B': 1.8, 'Al': 2.1, 'Ni': 2.0, 'Cu': 2.0, 'Se': 1.9, 'As': 2.05,
			'Ag': 2.1, 'Au': 2.1, 'Pd': 2.05, 'Pt': 2.05, 'Pb': 2.3, 'Hg': 2.05, 'others': 2.0}
	
	COLORS = {'N': 'blue',
			'C': 'gray',
			'O': 'red',
			'H': 'white',
			'Cl': 'yellow',
			'others': 'black'}
	
	BOND_COLOR = 'gray'
	
	def __init__(self, atom_radii=None, atom_colors=None):
		self.atom_radius = self.atom_radius = atom_radii if atom_radii else {k: self.ATOM_SCALE * v for k, v in self.Rvdw.items()}
		self.atom_color = self.atom_color = atom_colors if atom_colors else self.COLORS.copy()

class BallStick(style_config):
	pass

class Sticks(style_config):
	RADIUS = 0.2   
	BOND_RADIUS = RADIUS
	ATOM_RADII = RADIUS

	COLOR = 'gray'
	BOND_COLOR = COLOR
	ATOM_COLORS = COLOR
	
	def __init__(self):
		super().__init__(atom_radii={k: self.RADIUS for k in self.Rvdw.keys()}, atom_colors={k: self.COLOR for k in self.COLORS.keys()})

class SpaceFilling(style_config):
	ATOM_SCALE = 1.