"""
Mol2Mesh
Parser functions for chemical table files

Written by Ira Novianti
"""
bond_types = {1.0: 'single', 2.0: 'double', 3.0: 'triple'}

def open_gjf(filename):
	text = open(filename, "r").read()
	sections = text.split('\n\n')

	method = sections[0]
	title = sections[1]

	atoms = []
	atoms_str = sections[2].split('\n')[1:]
	for line in atoms_str:
		data = line.split(None)
		element = data[0]
		coor = [float(data[-3]),
				float(data[-2]),
				float(data[-1])]
		atoms.append({'element': element,
					'coor': coor})

	bonds = []
	bonds_str = sections[3].split('\n')
	for line in bonds_str:
		data = line.split(None)
		if len(data[1:]):
			atom1_idx = int(data[0]) - 1
			linked_atom = data[1:]
			for i in range(0, len(linked_atom), 2):
				atom2_idx = int(linked_atom[i]) - 1

				coor_1 = atoms[atom1_idx]['coor']
				coor_2 = atoms[atom2_idx]['coor']

				bonds.append({'coor_1': coor_1,
							  'coor_2': coor_2})
	return {'atoms': atoms, 'bonds': bonds}


def open_mol(filename):
	text = open(filename, "r").read()
	content = text.split("M  END")[0].split('\n')[3:]

	info = content[0].split(None, 3)
	n_atoms = int(info[0])
	n_bonds = int(info[1])

	atoms_bonds = content[1:-1]

	atoms = []
	atoms_str = atoms_bonds[:n_atoms]
	for line in atoms_str:
		data = line.split(None, 5)
		coor = [float(data[0]),
				float(data[1]),
				float(data[2])]
		element = data[3]
		atoms.append({'element': element,
					  'coor': coor})

	bonds = []
	bonds_str = atoms_bonds[-n_bonds:]
	for line in bonds_str:
		data = line.split(None, 3)
		atom1_idx = int(data[0]) - 1
		atom2_idx = int(data[1]) - 1

		bond_type = bond_types[float(data[2])]

		coor_1 = atoms[atom1_idx]['coor']
		coor_2 = atoms[atom2_idx]['coor']

		bonds.append({'coor_1': coor_1,
					  'coor_2': coor_2,
					  'bond_type': bond_type})

	return {'atoms': atoms, 'bonds': bonds}
