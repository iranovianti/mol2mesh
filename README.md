<img src="https://raw.githubusercontent.com/iranovianti/mol2mesh/main/figures/logo2.png" width="213" height="162">

-----------

# About
**Mol2Mesh** generates a [3D model](#3d-models) of any molecule from its [chemical table file](#chemical-table-files)

### Input (Chemical table file)
Example: [Dopamine.sdf](https://github.com/iranovianti/mol2mesh/blob/main/sample_files/Dopamine.sdf)
```
681
  -OEChem-05142402263D

 22 22  0     0  0  0  0  0  0999 V2000
   -2.2392    1.9626    0.0548 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3557   -0.5612    0.3868 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4081    0.2624    0.3445 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.1628   -0.0212   -0.6613 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7040   -0.1603   -0.3850 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9862    0.1008    0.6289 C   0  0  0  0  0  0  0  0  0  0  0  0
```

### Output (3D model)
  <img src="https://github.com/iranovianti/mol2mesh/blob/main/figures/Dopamine_stl_3Dviewer.gif" width="400" height="400"> <img src="https://github.com/iranovianti/mol2mesh/blob/main/figures/Dopamine_glb_3Dviewer.gif" width="400" height="400">

  3D models of dopamine `.stl` (left) and `.glb` (right) viewed in 3D viewer

# Basic usage
<a href="https://colab.research.google.com/drive/1wXtnuMb0TBxx6yrmy46tHzWMkTw6Mk3w"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"></a>
``` python
from mol2mesh import Mol2Mesh

#open chemical table file (.sdf, .mol, or .gjf)
input_file = 'sample_files/Dopamine.sdf' #path to file
molecule = Mol2Mesh(input_file)

#check the resulting 3D model of the molecule
molecule.show()

#save 3D model
molecule.save_mesh('Dopamine.stl')
```

### Parameters
``` python
Mol2Mesh(file_path, style='BallStick', res_a=25, res_b=15, multicov=False, name=None)
```
* `file_path`: Path to the chemical table file (sdf, mol, or gjf).
* `style`: Style configuration defining the radii and colors of atoms and bonds. It can be a [configuration class](https://github.com/iranovianti/mol2mesh/blob/main/style.py) or a string from the defined styles ['BallStick', 'Sticks', 'SpaceFilling'].

  <img src="https://raw.githubusercontent.com/iranovianti/mol2mesh/main/figures/styles.png" height="200" width="675">

* `res_a` and `res_b`: Resolution of an atom sphere and bond cylinder, i.e., the number of points defining each circle forming the sphere or cylinder. Increase the value for a smoother surface and decrease for a smaller file size.

  <img src="https://raw.githubusercontent.com/iranovianti/mol2mesh/main/figures/res_a.jpg" height="300" width="657">

* `multicov`: Visibility of multiple covalent bonds
  
  <img src="https://raw.githubusercontent.com/iranovianti/mol2mesh/main/figures/multicov.jpg" height="218" width="369">
  
* `name`: The name of the molecule to be stored in the resulting 3D model `.stl` file.

### Methods
* `show()`: Visualizes the molecule in 3D using Plotly.
* `save_mesh()`: Saves the molecule mesh as an `.stl` file. The resulting file can be used for 3D printing. `.stl` files do not store color.
* `get_mesh()`: Returns a combined array of meshes as np.array
* `show_mesh()`: Displays molecule meshes (vertices and edges) as wireframes using Plotly.
* `save_glb()`: Saves the molecule as a 3D model with colors in `.glb` format. Relies on trimesh library.

## Chemical table files
[Chemical table file](https://en.wikipedia.org/wiki/Chemical_table_file) (CT file) is a text file that describes a molecule. A CT file can store many information, but at the very least it should contain xyz coordinates of atoms that make up a molecule. There are many different formats, but the most common one is molfile which can have `.mol` or `.sdf` file extension.
### Getting chemical table files for known compounds
* [PubChem](https://pubchem.ncbi.nlm.nih.gov/): Most known small-molecule compounds are registered in PubChem. To get the sdf file, search for the molecule by name, and on its molecule page, go to `3D conformer` section to download the sdf file.

  <img src="https://raw.githubusercontent.com/iranovianti/mol2mesh/main/figures/get_sdf_pubchem.jpg" height="250" width="292">
* Wikipedia: If you don't know the exact name of the molecule, you can google the molecule and check whether it has a wikipedia page.  In the identifier section (usually located beneath the molecule's structure on the right side of the page), you'll find the PubChem CID, which will redirect you to the molecule's page on PubChem.
  
  <img src="https://raw.githubusercontent.com/iranovianti/mol2mesh/main/figures/get_sdf_wikipedia.jpg" height="300" width="200">

### Creating chemical table files
For custom molecules, you can use any 3D molecular visualization software. The most simple and acessible one is [MolView](https://molview.org) which allows you to download the `.mol` file of the molecule, however it doesn't offer energy minimization, so the 3D conformer for unknown compounds might look strange. Alternatively, you can use well-known software such as Avogadro, Chem3D, MarvinView, etc., and save the molecule as `.mol`, `.sdf`, or `.gjf`.

## 3D models
### Types
* `.stl` file is a human-readable file  of a 3D mesh that contains triangles (xyz coordinates) and their facet normals.
* `.glb` file is a representation of a 3D model written by GL Transmission Format (glTF).

### Usage
The resulting file can be used for educative or personal purposes, such as making a real 3D model by 3D printing or just for molecule visualization. Additionally, both `.stl` and `.glb` files can be imported into PowerPoint slides.
  
  <img src="https://raw.githubusercontent.com/iranovianti/mol2mesh/main/figures/Dopamine_powerpoint.gif">
