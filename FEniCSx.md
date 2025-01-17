# notes on FEniCSx

## Installation
- Currently developed actively, and api also changes
- 
- has complex mode and real mode; PDE solver library PETSc has either real or complex build, and FEniCSx has to be switched to the corresponding mode. Below command to check the available version and build, and install PETSc(Corresponding FEniCSx will be automatically installed). It is an idea to install both in different environment to allow switching.


'''
conda search -c conda-forge PETSc
conda install -c conda-forge petsc=3.21.6=complex_h9a830cf_0

'''

# Coding note
## When modifying the code written for real mode to run on complex mode
	- Wrap constant with PETSc.ScalarType [Ref](https://jsdokken.com/dolfinx-tutorial/chapter1/complex_mode.html#variational-problem/)
	```python
	rho = fem.Constant(domain, PETSc.ScalarType(2700.0))
	```
		- Expected error [ incompatible function arguments. The following argument types are supported: ]
	- conjugate tensor product either by ufl.conj or ufl.inner(second operand conjugated)	
		- Expected error [ ArityMismatch: Failure to conjugate test function in complex Form ]
		
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


## When modifying from one version to another
- ERROR [compute_integration_domains(): incompatible function arguments.]
  - Some versions needs explicit type definition for meshtags
```python
facet_tags = meshtags(mesh, mesh.topology.dim - 1,
		np.concatenate([facets_left]),
		np.array([0] * len(facets_left)).astype('int32'))
```




# Zakki
- [Jupyter lab Text Editor enhancement](https://github.com/krassowski/jupyterlab-lsp#highlight-references>):[Installation](https://jupyterlab-lsp.readthedocs.io/en/latest/Installation.html)

- Install relevant traum libraries and jupyter labextension using conda(used pip to install viskex at first, ending up with those traum libraries installed by pip and not recognized)[Tram for Jupyter Lab](https://www.kitware.com/trame-in-jupyterlab-a-unified-approach-for-web-apps-in-scientific-computing/)
- [How to set up traum for pyvista](https://tutorial.pyvista.org/tutorial/00_jupyter/index.html)
- [Different ways to slice in pyvista](https://docs.pyvista.org/examples/01-filter/slicing.html)
- [Difference between Complex mode and Real modeheawu](https://jsdokken.com/dolfinx-tutorial/chapter1/complex_mode.html)
- [FEniCSx tutorial for Sorbonne University](https://jsdokken.com/FEniCS23-tutorial/README.html)
- [FEniCSx official tutorial](https://jsdokken.com/dolfinx-tutorial/chapter1/complex_mode.html)



- [General modal analysis for fenics](https://fenics-solid-tutorial.readthedocs.io/en/latest/EigenvalueProblem/EigenvalueProblem.html)



- [COMSOL6.3-Shell and Plate introduction-](https://doc.comsol.com/6.3/doc/com.comsol.help.sme/IntroductionToStructuralMechanicsModule.pdf)
```
The Shell interface ( ) is intended for mechanical analysis of thin-walled structures. The formulation used in the Shell interface is a Mindlin–Reissner type, which means that transverse shear deformations are accounted for. It can be used for rather thick shells as well as thin ones. It is possible to prescribe an offset in a direction normal to a selected surface, which for example can be used when meshing imported geometries. The Shell interface also includes other features such as damping, thermal expansion, and initial stresses and strains. 
```
- [Structural Mechanics Theory(p.494 for elastic moduli)](https://doc.comsol.com/6.3/doc/com.comsol.help.sme/StructuralMechanicsModuleUsersGuide.pdf)

- [FEnicSx-Shell](https://fenics-shells.readthedocs.io/en/latest/demo/reissner-mindlin-clamped/demo_reissner-mindlin-clamped.py.html)


-[Jupytext to convert .py <-> .ipynb](https://jupytext.readthedocs.io/en/latest/using-cli.html): jupytext --to notebook notebook.py <-> jupytext --to py notebook.ipynb     
