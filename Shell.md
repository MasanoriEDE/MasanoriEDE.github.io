# Vector Analysis
## Vector Algebra
- Three types of vectors
    - Free vector: can be displaced parallel to itself
    - Sliding vecotr: can be displaced along its line of action
    - Bound vector: has the origin and destination
- a quantity characterized by a magnitude and a direction is not necessaryly a vector
    - **The quantity must obey the laws of vector algebra**
- Vector addition
    - commutative:``A+B=B+A``
    - associative: ``(A+B)+C=A+(B+C)=A+B+C``
- Vector subtraction
    - the same magnitude but with the opposite direction
- Projection of a vector onto an axis $u_0$
    - A projectino of a vector: ``A_u=A\cos\phi=A\cos(A,u_0)``
- Multiplication of a vector by a scalar
    - multiply the magnitude of the vector

## Bases and basis vectors
- a three-dimentional space can be defined by a set of three linearly independent vectors
    - oblique coordinate system(``(x^1,x^2,x^3)``) vs orthonormal coordinate(``(x_1,x_2,x_3)``)
- The position of a point M is uniquely determined by its radius vector relative to the fixed origin in the coordinate system
    - However, before making any calculations, we must introduce a coordinate sysytem to denote the coordinates of the point in the coorinate system and in a certain unit
- Orthonormal systems: a coordinate system whose basis vectors intersect at right angles 
- Curvilinear coordinate: a coordinate system whose  coordinate curves are not straint lines(<-> rectangular, oblique coordinates)
    - Polar coordinate is a curvilinear orthonormal coordinate except when elliptic

## Products of two vectors
- Scalar product(dot product): ``A\cdot B=|A||B|\cos(A,B)=A_1B_1+A_2B_2+A_3B_3``
    - commutative: ``A\cdot B=B\cdot A``
    - distributive: ``A\cdot (B+C)=A\cdot B+A\cdot C``
    - ``A\perp B \leftrightarrow A\cdot B=0``
- Vector product(cross product)
    - magnitude: ``|A||B|\sin(A,B)``
    - orientation: perpendicular to the plane of A and B, direction from which rotation from A->B appears to be counter-clockwise
    - non-commutative: ``A\times B=-B\times A``
    - ``A||B \leftrightarrow A\times B=0``
    - ``C=A\times B=\begin{vmatrix}
i_1 & i_2 & i_3\\
A_1 & A_2 & A_3\\
B_1 & B_2 & B_3\\
\end{vmatrix}
``
- Scalar triple product: ``V=(A\times B)\cdot C=|A\times B|C_{A\times B}=|A\times B|h``
    - the volume of the parallelepiped taken with the sign depending on whether the angle between ``C`` and ``A\times B`` is acute or obtuse
    - ``(A\times B)\cdot C=\begin{vmatrix}
i_1 & i_2 & i_3\\
A_1 & A_2 & A_3\\
B_1 & B_2 & B_3\\
\end{vmatrix}\cdot(C_1i_1+C_2i_2+C_3i_3)=
\begin{vmatrix}
C_1 & C_2 & C_3\\
A_1 & A_2 & A_3\\
B_1 & B_2 & B_3\\
\end{vmatrix}
``
    - commutative: ``(A\times B)\cdot C=(B\times C)\cdot A=(C\times A)\cdot B``
    - ``(A\times B)\cdot C=0`` if three vectors are coplanar(linearly dependent)
- Vector triple product: ``A\times (B\times C)``
    - perpendicular to A and lies on the B,C plane


## 4.3: Scalar Fields
### Level Surfaces
### Gradient and directional derivative
Gradient vector: ``grad \phi=i_i\frac{\partial \phi}{\partial x_1}+i_2\frac{\partial \phi}{\partial x_2}+i_3\frac{\partial \phi}{\partial x_3}=i_k\frac{\partial \phi}{\partial x_k}``
    - describes the inhomogeneity of the field
Directional derivative: ``\frac{d\phi}{dl}=l\cdot grad \phi`` where l is a unit vector

- Properties of the gradient
``\frac{d\phi}{dl}=|grad\phi|cos(l,grad\phi)``
    - The rate of change is the greatest in the direction of $grad \phi$: ``\bigg(\frac{d\phi}{dl}\bigg)_{max}=|grad\phi|=\sqrt{\bigg(\frac{\partial \phi}{\partial x_1}\bigg)^2+\bigg(\frac{\partial \phi}{\partial x_2}\bigg)^2+\bigg(\frac{\partial \phi}{\partial x_3}\bigg)^2}``
    - vecotr grad $\phi$ points at the direction that's normal to the level surface
- $\nabla$ ('del' or 'nabla') symbol is used 

## 4.4 Vector Fields
### Trajectory
### Flux
``A\cdot n dS=A_n dS`` or by integrating through the surface ``\int \int_SA\cdot n dS=\int\int_SA_ndS``
- let dS be an element of smooth surface S, Flux is the amount flowing through dS per unit time
    - note $n$ is the unit exterior normal to S
- dS defines an elementary tube of flow whose surface is formed by the trajectories of the vector field
- the flow amount can be calculated by the height of the tube``h=|v|dt\cos(v,n)=|v\cdot n|dt``
    - the flow amount is then: ``dQ=|v|cos(v,n)dS=v\cdot ndS``
    - by integrating over the whole surface: ``Q=\int \int _Sv\cdot n dS``
- if the fluid is incompressible, the mass m of the flow through the surface S per unit time is ``m=\rho Q``
- if the flux of a field of a volume V, enclosed by a surface S, is positive, then flow out > flow in: decrease in density
### Divergence
``\mathrm{div}A=\lim_{V\rightarrow0}\frac{1}{V}\int\int_SA\cdot ndS=\frac{1}{V}\int\int_V\bigg(\frac{\partial A_1}{\partial x_1}+\frac{\partial A_2}{\partial x_2}+\frac{\partial A_3}{\partial x_3}\bigg)dV=\frac{\partial A_1}{\partial x_1}+\frac{\partial A_2}{\partial x_2}+\frac{\partial A_3}{\partial x_3}``
- the flux of A through S divided by V is the average strength of the sources and sinks inside V
    - its limit as the V and S shrinks to a point M: divergence
    - the definition of div A is independnet of the choice of coordinate system
- field  div A does not exist for every field A,
        - but only exists at every point where the components and thier derivatives are continuous
- From Gauss' theorem, the surface integration may be carried over to the volume integration
    - this volume integration of components yields the definition of divergence of A in a rectangular coordinates
- Also, by rewriting Gauss' theorem with divergence, divergence theorem is obtained: ``\int\int\int_V\mathrm{div} A dV=\int\int_S A\cdot ndS``
- ``\mathrm{div}A`` is a sca
- lar product of $\nabla$ and $A$: ``\mathrm{div} A=\frac{\partial A_k}{\partial x_k}``


# Mode
Compliance=Displacement/force
mechanical impedance=1/mobility=force/velocity
accelarance=acceleration/force

# Shell Physics


free vector
bound vector

## Vector Fields "Vector and tensor analysis with applications"
Trajectories: A curve whose tangent at every point has the same direction as a vector field ``A=A(r)``
Flux: Let S be a two-sided piecewise-smooth surface in a vector field ``A(r)``, whose elment is called dS. ``n`` is a unit vector normal to dS. Flux is the ``A\cdot n dS=A_ndSf``
### Divergence
- Given any point M in a vector field A=A(r), let S be an arbitrary closed surface S surrounding M, and enclosing a volume V
    - ``\frac{1}{V}``
- the limit of V and S is called the divergence of the field A


## [4th order PDE: vibrating beam(Faculty of Kahn)](https://www.youtube.com/watch?v=eQ2MqMVc_3M&list=PLdgVBOaXkb9Ab7UM8sCfQWgdbzxkXTNVD&index=14)
Equation of vibration: ``\frac{\partial^2 u}{\partial t^2}=-\alpha^2\frac{\partial^4u}{\partial x^4}``
``\alpha-\frac{EI}{\mu}`` :flextual rigidity
- E:Young modulus, I:area moment of inertia w.r.t x-axis, $\mu$: mass per length
Boundary condition: ``u=0``
-> Bending moment of the beam:``\frac{\partial^2 u}{\partial x^2}=0`` (from structural mechanics point of view)

Initial conditions in the beam:
- $u(x,0)=u_0(x)$
- $\frac{\partial u }{\partial t}|_{t=0}=v_o(x)$

PDE: 4th order in x--4 initial conditions in x
    - $u(x,t)=u(L,t)=0$
    - ``\frac{\partial^2 u}{\partial x^2}|_{x=0}=\frac{\partial^2 u}{\partial x^2}|_{x=L}=0``
PDE: 2th order in t--2 initial conditions in t
    - ``u(x.0)=u_o(x), \frac{\partial u}{\partial t}|_{t=0}=v_0(x)``

- PDE and boundary conditions are homogeneous: use **separation of variables**--let $u(x,t)=X(x)T(t)$
    - then PDE is: $X\frac{d^2T}{d^2t}=-\alpha^2T\frac{d^4 X}{dx^4}$
    - or $\frac{1}{T}\frac{d^2 T}{dt^2}=\frac{-\alpha^2}{X}\frac{d^4X}{dx^4}$, which has to be constant
    - ![](assets/17400663908160.jpg)



## Formulation of equations of motion: "Introduciton to finite element vibration analysis"
Weak form for CPT:..?
``D\int_\Omega\nabla ^4wvd\Omega+\rho h\ddot w=f``
``D\int_\Omega\nabla^2w\nabla^2vd\Omega+D\bigg[\int_S \bigg\{ v(\nabla(\nabla w)n) - (\nabla w)(\nabla v n)\bigg\}\bigg] + \rho h\ddot w=f``



## Strong and weak forms
Strong forms: governs the system of equations, requires strong continuity on the dependent field variables--differentiable up to the order of the partial differential equations(=strong form)
- Obtaining exact solution for a strong form is usually very difficult
  - The finite difference method(FDM): anapproximated solution to strong form, but only works for simple and regular geometry and boundary conditions

Weak form: usually created using one of the following methods, a weaker continuity on the field variables, produces a set of discretized system equations that give much more stable and accurate results especially for complex geometry
- Energy principles: special form of variational principle, suited to mechanics of solids and structures
- Weighted residual methods: general mathematical tool for all kinds of PDE


  
``\kappa hG\nabla^2w+\kappa h G\Phi-
\rho h\frac{\partial ^2w}{\partial t^2}=-F``
OO
``\gamma_{bs}\big\{\big\}``


``\nabla^2_{\eta\xi}W+\bar\Phi_{\eta\xi}+\gamma_{bs}\Omega^2W=0``
``r_{o}^2\gamma_{bs}\bigg\{
\frac{1}{2(1+v)}\nabla^2_{\eta\xi}\Psi_x+
\frac{1}{2(1-v)}\frac{\partial\bar\Phi_{\eta\xi}}{\partial\eta}\bigg\}+
(r_o^2\gamma_{bs}\Omega^2-1)\Psi_x-
\frac{\partial W}{\partial \eta}
=0``
``r_{o}^2\gamma_{bs}\bigg\{
\frac{1}{2(1+v)}\nabla^2_{\eta\xi}\Psi_y+
\frac{1}{2(1-v)}\frac{\partial\bar\Phi_{\eta\xi}}{\partial\xi}\bigg\}+
(r_o^2\gamma_{bs}\Omega^2-1)\Psi_y-
\frac{\partial W}{\partial \xi}
=0``
where$$\bar\Psi_{\eta\xi}=\frac{\partial\Psi_x}{\partial\eta}+\frac{\partial\Psi_y}{\partial\xi}=\frac{1}{1+v}(\bar M_x+\bar M_y)$$


## 
![](assets/17398833766940.jpg)

![](assets/17398850913378.jpg)


## The Finite Element Method
stress
strain
deformation
displacement
Solid mechanics
Structural mechanics
Static<Dynamic
Elastic,Elasticity
Plastic,Plasticity
linear relationship
linear elastic
anisotropic<isotropic


deflection:w
rotation about x:theta_x
rotation about y:theta_y

hoge


## Solid Mechanics Theory
### Tensor notation

### Generalized Plane Strain(P.269)

## Vectors and tensors

$\epsilon$: Physical Euclidean space, R^3(set of three real numbers)
- O: origin of $\epsilon$

ξ1, ξ2 and ξ3 are also called the coordinates of M in the coordinate system defined by O and (⃗ı1,⃗ı2,⃗ı3). 






Invariant Quantities: qualities that remain unchanged regardless of how the coordinate system is set(scalars)
Covariant Quantities: qualities that follow the coordinate system:when the coordinate is stretched or rotated, the quality follows the same 
	- spacial derivatives(gradients etc.) depend on the system scale and rotation
	- forces depends on the orientation
Contravariant Qualities: qualities that oppose the change in the coordinate system
	- displacements and velocities shrinks when the coordinate is streched and also rotates against the coordinate

Contravariant basis:superscript
Covariant basis:subscript

Components of a vector in the coontravariant basis: covariant components?
Components of a vector in the covariant basis: contravariant components?


Metric tensor: dot product of covariant-covariant components, or contravariant-contravariant components
* Why "metric"--dot product calculates the length of vectors ie distances
* If covariant and contravariant are orthgonal, there's one metric tensor between them

$$g_{mn}=i_m\dot i_n$$

Metric can be used to exchange between contravariant and covariant
- the same goes for component-wise exchange via components of the metric tensor



*To distinguish surface tensors from 3D tensors, surface tensors are denoted with a number of underbars corresponding to their order*
first fundamental form of the surface: one under bar, metric tensor of the tangent plane, can be used to convert covariant components into contravariant ones, useful to express surface integrals,
second fundamental form of the surface: two under bar


Dual vector: can be described by covariant component, takes contravariant vector, returns a scalar


The total rank of V is m=p+q: (p,q)-tensor, where
	p:contravariant rank
	q:variant rank
Tensor operation such as Summation, Scalar Multiplication, Linear Combination --they do not change the rank of the tensorst

Cordinates

Tensor: an object that generalizes vector to a higher dimension, that's invariant under a change of coordinate systems, with components that change according to a special set of mathematical formulae 
	
Tensor product:


Tangent space:A vector space that tangentially touches a smooth surface M at point T, described as $T_{p}M$,represents physical directions of movement
Dual space:The dual space of the tangent space at p($T^{*}_{p}M$) is the set of all linear functionsthat take all tangent vectors v and return a scalar,$\omega:T_{p}M->R$
Metric tensor is a bilinear map.represents the measurements applied to vectors

Tangent Space: The space of vectors that describe directions at a point.
Dual Space: The space of covectors (linear maps) that "act on" tangent vectors to produce scalars.
Metric Tensor: The tool that connects these spaces, allowing transformations between covariant and contravariant components.

Injective mapping: a function that maps elements from one space to another in a way that preserves uniqueness--no two discinct elements in the original space are mapped to the same element in  the target space.

The metric tensor $g_{ij}$ 


Chart:A smooth injective mapping from one bounded open subset of R, $Omega$ into another $\epsilon$

Invariant quantity: a quantity that does not depend on a particular choice of coordinate system. Velocity doesn't depend on the coordinate, but it's components do depend on the coordinate.







![](assets/17407606966984.jpg)


![](assets/17410846879231.jpg)
(ref)[http://inspire.starfree.jp/analyticalMechanics/introduction.html]

https://physnotes.jp/math/par_tot/