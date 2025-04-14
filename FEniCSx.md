# notes on FEniCSx

## Installation
- Currently developed actively, and api also changes
- 
- has complex mode and real mode; PDE solver library PETSc has either real or complex build, and FEniCSx has to be switched to the corresponding mode. Below command to check the available version and build, and install PETSc(Corresponding FEniCSx will be automatically installed). It is an idea to install both in different environment to allow switching.


'''
conda install -c conda-forge fenics-dolfinx==0.9.0 mpich pyvista
conda search -c conda-forge PETSc
conda install -c conda-forge petsc=3.22.3=complex_h1434fe3_0
(for 0.8.0, conda install -c conda-forge petsc=3.21.6=complex_h9a830cf_0)
conda install -c conda-forge meshio
conda install -c conda-forge jupyterlab h5py jupyterlab_widgets
conda install -c conda-forge matplotlib ipywidgets plotly
cd viskex
python3 -m pip install '.[tutorials]' 
source /opt/homebrew/Caskroom/miniconda/base/envs/FEMBEM25Mar 
which python
pip install bempp-cl
'''

# Coding note

## tabulation
evaluate the basis function at a certain set of points in the reference element by using tabulation.
The first argument: the number of spatial derivatives of the basis functions we want to compute
The second argument: a set of input points to evaluate the basis functions at as a numpy array of shape (num_points, reference_cell_dimension)
```python
points = np.array([[0.0, 0.1], [0.3, 0.2]])
values = element.tabulate(0, points)
```
Outputs the array of shape (num_spatial_derivatives+1, num_points, num_basis_functions)

In some cases, one might want to use a lower accuracy for tabulation of basis functions to speed up computations. This can be changed in basix by adding dtype=np.float32 to the element constructor.
Basix also allows other extra options to tweak finite elements--try changing the options and see the shape of the shape functions


## Mapping
As we have already seen, we can describe any cell in our subdivided domain with a mapping from the reference element. However, as we want to integrate over each element individually, we need to map the basis functions to and from the reference element. We call this map: 
.
```
What property does the covariant Piola mapping preserve?
We observe that the covariant Piola map maps normals to normals, since a 0 tangential component is mapped to a 0 tangent.
```

## Entity map
[This](https://fenicsproject.org/blog/v0.9.0/#mixed-assembly) guides you how to create entity map to assemble across subdomain and main domain.
Entity map is the reverse dictionary(sub->parent domain) of the parent->sub mapping obtained as the 2nd return of mesh.create_submesh.


## creating entities and connectivity/ obtaining boundary facets
```python
mesh.topology.create_entities(mesh.topology.dim - 1)
gmshmesh.topology.create_connectivity(mesh.topology.dim - 1, mesh.topology.dim)
boundary_facets = dolfinx.mesh.exterior_facet_indices(mesh.topology)
```

```
When a mesh is created in DOLFINx we only compute the unique global numbering and ownership of cells and vertices. However, one can compute this numbering for any sub-entity by calling mesh.topology.create_entities(j). For instance computing the ownership and numbering of the facets can be done with mesh.topology.create_entities(mesh.topology.dim-1).
```

```
When a mesh is created in DOLFINx, the relation between the cells and the vertices is computed. However, one can compute the relationship between any set of entities in the mesh with dolfinx.mesh.create_connectivity(i, j) where i is the dimension of the index to map from and j is the dimension of the entities to map from. As an example, calling mesh.topology.create_connectivity(mesh.topology.dim-1, mesh.topology.dim) will compute the relationship between all facets in the mesh and their cells. The inverse map (cell-to-facet) is computed with mesh.topology.create_connectivity(mesh.topology.dim, mesh.topology.dim)-1
```

With these two preliminaries computed, DOLFINx can determine which facets that are connected to a single cell only.

## 
Another feature of variational formulations is that the test function v is required to vanish on the parts of the boundary where the solution u is known. See for instance [LM19].



This simple but very powerful method for constructing test problems is called the method of manufactured solutions. First pick a simple expression for the exact solution, plug into the equation to obtain the right-hand side (source term 
f
). Then solve the equation with this right hand side, and using the exact solution as boundary condition. Finally, we create a program that tries to reproduce the exact solution.



Entities vs degrees-of-freedom
Entities: A mesh entity in FEniCS is either a vertex, an edge, a face, or a cell (triangle or tetrahedron)
Degrees-of-freedom: mesh nodes and interpolation points
ct vs ft
ct: volumetric element eg)hexahedron
ft: surface element eg)quadrilateral

## Creating a source on a pointO
```
Next, we create a function containing the exact solution (which will also be used in the Dirichlet boundary condition) and the corresponding source function for the right hand side. Note that we use ufl.SpatialCoordinate to define the exact solution, which in turn is interpolated into uD and used to create the source function
```

## Defining a boundary condition for mixed elements
V.sub(0).sub(0) allows accessing the item of a list that's a part of the mixed element.
```python
    clamped_dofsW = fem.locate_dofs_topological(V.sub(0), facet_dim, clamped_facets)
    bcs = [fem.dirichletbc(PETSc.ScalarType(0.0), clamped_dofsW, V.sub(0))]
    clamped_dofsT = fem.locate_dofs_topological(V.sub(1), facet_dim, clamped_facets)
    bcs += [fem.dirichletbc(PETSc.ScalarType(0.0), clamped_dofsT, V.sub(1).sub(0))]
    bcs += [fem.dirichletbc(PETSc.ScalarType(0.0), clamped_dofsT, V.sub(1).sub(1))]
```

## Avoid mixed elements for coupled PDEs
When assembling the matrix blocks, mixed elements block is created to construct its own diagonals, ending up with a different order of dofs--BCs are applied to wrong dofs in the end.


## When porting codes (from FEniCSx version, from different python environment)
### When modifying the code written for real mode to run on complex mode
	- Wrap constant with PETSc.ScalarType [Ref](https://jsdokken.com/dolfinx-tutorial/chapter1/complex_mode.html#variational-problem/)
	
```python
	rho = fem.Constant(domain, PETSc.ScalarType(2700.0))
```		

	- Expected error [ incompatible function arguments. The following argument types are supported: ]
	- conjugate tensor product either by ufl.conj or ufl.inner(second operand conjugated)	
		- Expected error [ ArityMismatch: Failure to conjugate test function in complex Form ]
	- [# Similarly, if we want to use the function `ufl.derivative` to take derivatives of functionals, we need to take some special care. As `ufl.derivative` inserts a `ufl.TestFunction` to represent the variation, we need to take the conjugate of this to be able to use it to assemble vectors.](https://github.com/jorgensd/dolfinx-tutorial/blob/main/chapter1/complex_mode.py)
		
		
	- Defining boundary condition with constants, below failed in solving the problem due to an incompatible type error
	
```python
bcs = [
    fem.dirichletbc(np.zeros((2,)), left_dofs, V),
    fem.dirichletbc(np.zeros((2,)), right_dofs, V),
]
```
	- Changing to numpy to PETSc.ScalarType resolved the issue
	
```python
bcs = [
    fem.dirichletbc(PETSc.ScalarType([0, 0]), left_dofs, V),
    fem.dirichletbc(PETSc.ScalarType([0, 0]), right_dofs, V),
]
```
PETSc.ScalarType wraps the constant source on the right hand side. 
This is because we want the integration kernels to assemble into the correct floating type.

- When plotting, some values have to be set back to real values before plot

```python
grid.point_data["Deflection"] = u.x.array.reshape(-1, 3).real
```


### When modifying from one version to another
- ERROR [compute_integration_domains(): incompatible function arguments.]
  - Some versions needs explicit type definition for meshtags
```python
facet_tags = meshtags(mesh, mesh.topology.dim - 1,
np.concatenate([facets_left]),
		np.array([0] * len(facets_left)).astype('int32'))
```

[FEniCSx API changne:Misc](
https://github.com/jorgensd/dolfinx-tutorial/pull/135/commits/af935961f2dab5b7d6fd7fff5c992f33adf7aca4#diff-9efba6744d3fa48d088cf6c111c2ea3071b18b615ab9140aa7065cd40ea57f30)

[FEniCSx API change:Function](https://github.com/FEniCS/dolfinx/pull/2818)

### When getting the value of mesh
```python
display(kappa.x.array[:])
```

### When importing mesh to FEniCSx
- gmshio seems to fail to read hexahedron...Unknown cell type-KeyError: np.int32(6)
  - works fine with tetrahedron
- cell:Volume physical group number
- facet_tags: Surface physical number


- viskex can only plot dolfinx.mesh.Mesh object, not meshio mesh object
  - so better to write in xdmf and import from dolfinx <- not necessaryly cf)dolfinx.io.gmshio.model_to_mesh
  - also, this process helps separating different elemtns of the mesh, such as hexa, triangle, quad and line
  - Importantly, xdmf can not handle mixed elements, so different elements have to be stored separately

### Setting boundary condition
```python
def Omega_0(x):
    return x[1] <= 0.5


def Omega_1(x):
    return x[1] >= 0.5
```
Note that both functions use a ≤ or ≥, as FEniCSx will evaluate each cell at all of the vertices, and thus has to return True for all vertices aligned with the interface to be marked properly.

### Defining material property
- when density is a part of the sound source term, the entire volume domain will be integrated without density in the integrand--move the denisty to the LHS so that density can be defined as a part of volume integration where different density can be defined on different part of the volume.

## Dealing with Mesh
### Check the shape of the imported mesh
```python
msh, cell, facet_tags=gmshio.read_from_msh('./data/MembraneAbsorber.msh',MPI.COMM_WORLD,0,gdim=3)
topology, cell_types, geometry = plot.vtk_mesh(msh, 3)
grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
p = pyvista.Plotter()
p.add_mesh(grid, show_edges=True)
p.view_xy()
p.show_axes()
p.show()
```
### Check the size of the imported mesh
create measures for surface and volume, and then integrate over the region to calculate the area and volume--making sure the size of the mesu is correct
```python
msh, cell, facet_tags=gmshio.read_from_msh('./data/MembraneAbsorber.msh',MPI.COMM_WORLD,0,gdim=3)
display(cell.values,facet_tags.values)
ds=Measure('ds',domain=msh,subdomain_data=facet_tags)
dx=Measure('dx',domain=msh,subdomain_data=cell)
AreaTop=fem.assemble_scalar(fem.form(1.0*ds(2)))
VPoro=fem.assemble_scalar(fem.form(1.0*dx(1)))
VAir=fem.assemble_scalar(fem.form(1.0*dx(100)))

display('Expecting:0.8x0.8=0.064:',AreaTop)
display('Expecting:0.05x0.8:0.04',VPoro)
display('Expecting:0.75x0.8:0.6',VAir)
```

## Solving a coupled PDEs
### Monolithic Direct Solver
- Assemble all the unknowns in a single large matrix to be solved simulatneously using a direct solver such as LU decomposition
- simple to implement, free from convergence issues
- Expensive for large problems, memory intensive

### Nested (block) matrix solver
- Split the matrix into blocks and solve each separately and iteratively
- One unknown is solved while other unknowns are treated as preconditioners?
- Memory efficient, more suitable for large problems, specialized solvers can be employed for each block
- Complex implementation, could be slower for small problems

### Iterative solver (Schur complement)
- Instead of solving the entire system directly, solves for one unknown at a time using Schur complement
- Highly efficient and scalable for larg systems
- Complex implementation and the performance is reliant on the choice of preconditioners

### Non-block direct solver
- Solves like monolithic block direct solver, but the structure of the system is not preserved and is treated as a single dense system
- easy to implement, as the whole matrix is solved with a single solver
- less efficient for large problem

# Zakki

## Jupyter
- [Jupyter lab Text Editor enhancement](https://github.com/krassowski/jupyterlab-lsp#highlight-references>):[Installation](https://jupyterlab-lsp.readthedocs.io/en/latest/Installation.html)

- Install relevant traum libraries and jupyter labextension using conda(used pip to install viskex at first, ending up with those traum libraries installed by pip and not recognized)[Tram for Jupyter Lab](https://www.kitware.com/trame-in-jupyterlab-a-unified-approach-for-web-apps-in-scientific-computing/)
- [How to set up traum for pyvista](https://tutorial.pyvista.org/tutorial/00_jupyter/index.html)
- [Different ways to slice in pyvista](https://docs.pyvista.org/examples/01-filter/slicing.html)
- [Jupytext to convert .py <-> .ipynb](https://jupytext.readthedocs.io/en/latest/using-cli.html): jupytext --to notebook notebook.py <-> jupytext --to py notebook.ipynb     
- [IBM:Jupyter Markdown Memo](https://www.ibm.com/docs/en/watson-studio-local/1.2.3?topic=notebooks-markdown-jupyter-cheatsheet)
- [Code formatter](https://jupyterlab-code-formatter.readthedocs.io/configuration.html)

## COMSOL
- [COMSOL6.3-Shell and Plate introduction-](https://doc.comsol.com/6.3/doc/com.comsol.help.sme/IntroductionToStructuralMechanicsModule.pdf)
```
The Shell interface ( ) is intended for mechanical analysis of thin-walled structures. The formulation used in the Shell interface is a Mindlin–Reissner type, which means that transverse shear deformations are accounted for. It can be used for rather thick shells as well as thin ones. It is possible to prescribe an offset in a direction normal to a selected surface, which for example can be used when meshing imported geometries. The Shell interface also includes other features such as damping, thermal expansion, and initial stresses and strains. 
```
- [Structural Mechanics Theory(p.494 for elastic moduli)](https://doc.comsol.com/6.3/doc/com.comsol.help.sme/StructuralMechanicsModuleUsersGuide.pdf)

- [Equation View](https://www.comsol.com/support/learning-center/article/Viewing-and-Accessing-the-Equations-and-Variables-for-Physics-Feature-Nodes-30761/122)

## FEniCSx
- [FEniCSX workshop: First Step](https://jsdokken.com/FEniCS-workshop/src/deep_dive/expressions.html)
- [FEniCS wokshop in India](https://github.com/iitrabhi/fenics-workshop/tree/main)
- [Difference between Complex mode and Real mode](https://jsdokken.com/dolfinx-tutorial/chapter1/complex_mode.html)
- [FEniCSx 2022:Helmholtz](https://jsdokken.com/fenics22-tutorial/helmholtz.html)
- [FEniCSx tutorial for Sorbonne University](https://jsdokken.com/FEniCS23-tutorial/README.html)
- [FEniCS workshop tutorial](https://jsdokken.com/FEniCS-workshop/src/unified_form_language/ufl_forms.html)
- [FEniCSx official tutorial](https://jsdokken.com/dolfinx-tutorial/chapter1/complex_mode.html)
- [Jorgen's Gist](https://gist.github.com/jorgensd)
- [General modal analysis for fenics](https://fenics-solid-tutorial.readthedocs.io/en/latest/EigenvalueProblem/EigenvalueProblem.html)
- [FEnicSx-Shell](https://fenics-shells.readthedocs.io/en/latest/demo/reissner-mindlin-clamped/demo_reissner-mindlin-clamped.py.html)
- [DOLFINx document](https://docs.fenicsproject.org/dolfinx/v0.9.0/python/demos/demo_biharmonic.html)
- [FEniCS shell(before FEniCSx)](https://fenics-shells.readthedocs.io/en/latest/demo/kirchhoff-love-clamped/demo_kirchhoff-love-clamped.py.html)
- [Defining subdomain after importing gmsh](https://jsdokken.com/dolfinx-tutorial/chapter3/subdomains.html)
- [FeniCSx Handson Examples](https://fenics-handson.readthedocs.io/en/latest/helmholtz/doc.html)
- [Electromagnetic demos](https://mikics.github.io/gsoc-jupyterbook/chapter1/demo_scattering_boundary_conditions.html)
- [Electromagnetic: PML and scattering](https://mikics.github.io/)
- [Example by UACh](https://uwemuehlich.github.io/FEniCSxDocumentedExamplesUACh/html/)
- [UKAN Helmholtz Equation](https://mhamilt.github.io/ca-knowledgebase.github.io/tutorials/fem/fenicsx-helmholtz-tutorial.html)
- [FEniCSX time harmonic wave equations](https://github.com/samuelpgroth/waves-fenicsx)
- [Stefano FEniCSx Tutorial: Helmholtz](https://computational-acoustics.gitlab.io/website/posts/31-intro-to-fenics-part-2/)
- [SciFEM](https://scientificcomputing.github.io/scifem/examples/point_source.html)
- [FEniCSx-pctools: FEnICSx preconditioning tools](https://gitlab.com/rafinex-external-rifle/fenicsx-pctools)

- [Newfrac](https://newfrac.github.io/fenicsx-fracture/notebooks/linear-elasticity/00-Mesh.html)
- [Elasticity training at Sorbonne University](https://newfrac.gitlab.io/newfrac-fenicsx-training/02-finite-elasticity/finite-elasticity-I.html)
- [ADIOS4FEniCSx](https://jsdokken.com/adios4dolfinx/docs/meshtags.html)
- [Solid Mechanics coupled](https://solidmechanicscoupledtheories.github.io/intro.html)
- [Mixed domain demos](https://github.com/jpdean/mixed_domain_demos)
- [FEniCSx based couped hydrogen transport-heat transfer simulation](https://festim.readthedocs.io/en/fenicsx/index.html)
- [Discontinuity at interface using mixed domains](https://fenicsproject.discourse.group/t/discontinuity-at-interface-using-mixed-domains/15040)
- [External operator for expressing non analyticals without ufl](https://github.com/a-latyshev/dolfinx-external-operator)
- [Random lecture notes](https://a654cc05c43271a5d22f-f8befe5e0dcd44ae0dccf352c00b4664.ssl.cf5.rackcdn.com/course/lectures/)
- [3D Thermo-poroelastplastic numerical simulations](https://github.com/Matt-L-McLean/poromechanics)
- [Hadrien p-FEM](https://github.com/HadrienLMS/Helmholtz_1D_pFEM/blob/main/Helmholtz_1D_pFEM.ipynb)
- [FEniCSx-preCICE adapter](https://github.com/precice/fenicsx-adapter)

- [Gmsh official document](https://gmsh.info/doc/texinfo/gmsh.html)
- [Gmsh with python reference with different topology](https://bbanerjee.github.io/ParSim/fem/meshing/gmsh/gmsh-meshing-for-code-aster/)
- [Gmsh with python reference](https://bbanerjee.github.io/ParSim/fem/meshing/gmsh/quadrlateral-meshing-with-gmsh/)
-[Gmsh Python API basics blog](https://neph.altervista.org/tutorial-gmsh-python-api-basics-mesh-creation/?doing_wp_cron=1734198489.9841949939727783203125)
- [Gmsh Python API for FEnics](https://jsdokken.com/src/tutorial_gmsh.html)
- [Gmsh on grasshopper](https://discourse.mcneel.com/t/gmsh-in-grasshopper/173125/8)
- [Johns Hopkins University Nguyen lab](https://me.jhu.edu/nguyenlab/doku.php?id=fenicsx)
- [cardiac mechanics with lumped model](https://github.com/marchirschvogel/ambit)
- [Dolfiny: FEniCSx high level wrapper](https://github.com/fenics-dolfiny/dolfiny)
- [Viskex tutorials](https://viskex.github.io/tutorials/02_garda/tutorial_garda_dolfinx.html)

- [PETSc](https://petsc.org/release/)

-[DefElement:an encyclopedia for finite element definitions](https://defelement.org/ciarlet.html)

## Bempp
- [DML modelling](https://www.diyaudio.com/community/threads/application-of-impulse-excitation-for-dml-design-and-analysis.383567/page-16)

## Physics
- [COMSOL CYCLOPEDIA: Eigenfrequency Analysis](https://www.comsol.com/multiphysics/eigenfrequency-analysis)

## Misc
- [Jurgen Bathe FEM lecture series](https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://www.youtube.com/watch%3Fv%3DoNqSzzycRhw&ved=2ahUKEwiIxZ-o8YmMAxWBUUEAHUwLDJQQtwJ6BAgREAI&usg=AOvVaw27Ma-zU3FPyxlYxBrOYHRN)
- [deal.II FEM](https://www.dealii.org/current/doxygen/deal.II/step_60.html)
- [FFmpeg](https://trac.ffmpeg.org/wiki/AudioChannelManipulation)

```math
V_{sphere} = \frac{4}{3}\pi r^3
```



# Everything else
## Studio design

## Acoustics
- [Siemens Modal Analysis](https://community.sw.siemens.com/s/article/dynamic-stiffness-compliance-mobility-and-more)
- [Siemens FRF](https://community.sw.siemens.com/s/article/what-is-a-frequency-response-function-frf)
- [Siemens Autopower function](https://community.sw.siemens.com/s/article/the-autopower-function-demystified)
- [Siemens Natural frequency](https://community.sw.siemens.com/s/article/Natural-Frequency-and-Resonance)
- [GA web model](https://interactiveacoustics.info/html/GA_IS_isRectangular.html)
- [Wave synthesis](https://www.osar.fr/notes/waveguides/)

## Manufactures
- [GRAS:microphone](https://www.grasacoustics.com/products/measurement-microphone-sets)


## Composers/Sound designers
- [Brian Oliveira](https://www.briandoliveira.com/case-studies)

## Rhino
- [Setting up Camera view](https://blinksandbuttons.net/how-to-set-camera-view-in-rhino/)

## Web template
- [1](https://lazyren.github.io/)
- [2](https://qwtel.com/resume/)
- [3](https://daniel.arneam.com/)
- [4](https://tseknet.com/about)
- [5](https://shawnyeager.com/)

## Shelves
- [gutenberg](https://www.gutenberg.org/)
