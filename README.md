Nonlinear Thermal Conduction {#nonlinear-thermal-conduction .unnumbered}
============================

This project is about practical aspects of solving sparse systems of
nonlinear equations $$\label{nonlin}
\mathbf{F}(\boldsymbol{U}) \ = \ \boldsymbol{0}$$ using the inexact
Newton method. In particular we will use parallel Newton--CG and look at
a case study on nonlinear thermal combustion in a self--heating medium
in a partially insulated square domain. This is the same model as the
one considered in the first project. This problem arises when
investigating critical parameter values for underground repositories of
self--heating waste which are partially covered by buildings. Beyond
certain critical parameter values the steady state solutions may become
very large or even unbounded leading to an explosion in the medium. For
more information see Greenway & Spence [@GrSp:85] and Adler [@Ad:83].

Model formulation {#model-formulation .unnumbered}
-----------------

We are concerned with finding solutions of the (dimensionless) nonlinear
thermal conduction equation $$\label{pde}
\mathcal{F}(u)[\lambda,\beta] \ := \ - \ \frac{\partial^2 u}{\partial x_1^2} \ - \
\frac{\partial^2 u}{\partial x_2^2} \ - \ g(u)[\lambda,\beta] \ = \ 0
\qquad \text{in } \Omega := [0,1]^2,$$ subject to the following mixed
boundary conditions on $\partial
\Omega$:

$$% \label{bc}
\hspace*{-2cm}
\begin{array}{c}
\text{Dirichlet} \left\{\begin{array}{rcll}
u & = & 0 & \quad \text{if } \ x_1 > \delta \ \text{ and } \ x_2 = 1
\end{array}\right.\\[1ex]
\text{Neumann} \left\{\begin{array}{rcll}
\displaystyle \frac{\partial u}{\partial x_1} & = & 0 & \quad \text{if } \ x_1
= 0 \ \text{ or } \ x_1 = 1\\[2ex]
\displaystyle \frac{\partial u}{\partial x_2} & = & 0 & \quad \text{otherwise}
\end{array} \right.
\end{array}$$

for $0 \le \delta < 1$. Note that these (more realistic) boundary
conditions differ from those from Assignment 1. The nonlinear functional
$$\label{arrhen}
g(u)[\lambda,\beta] \ := \ \lambda \; \exp\left( \frac{u}{1+\beta u} \right)$$
is the Arrhenius reaction rate which depends on the two parameters
$\lambda$ and $\beta$.

The relationship between the dimensionless variables in
([\[pde\]](#pde){reference-type="ref" reference="pde"}),
([\[arrhen\]](#arrhen){reference-type="ref" reference="arrhen"}) and the
physical quantities is given by Adler [@Ad:83]. We note that $u$ is a
dimensionless temperature excess and $\lambda$ the Frank--Kamenetskii
parameter; $\beta$ is the dimensionless activation energy and $\delta$
the dimensionless half--width of an insulating strip.

Finite Difference Discretisation {#finite-difference-discretisation .unnumbered}
--------------------------------

As for the linear Poisson equation, suppose that $\Omega$ is subdivided
into small equal squares of size $h \times h$ with $h := 1/m$. The
interior and boundary nodes of the mesh are given by
$\boldsymbol{x}_{i,j} := (ih,jh)$, with $i,j =
0,\ldots,m$. Let $m_\delta := [\delta m] \in \mathbb{N}$ be the integer
part of $\delta m$ (i.e. largest integer $m_{\delta}\le \delta m$). Then
$m_\delta h
\le \delta$ and $(m_\delta + 1)h > \delta$, and so all the points
$\boldsymbol{x}_{i,m}$ with $i > m_\delta$ lie on the Dirichlet
boundary.

The finite difference approach to the above problem is now to seek
approximations $u_{i,j}$ to $u(\boldsymbol{x}_{i,j})$ such that
$$- \frac{u_{i-1,j} - 2 u_{i,j} + u_{i+1,j}}{h^2} \ - \
\frac{u_{i,j-1} - 2 u_{i,j} + u_{i,j+1}}{h^2} \ - \
g(u_{i,j}) \ = \ 0.$$ We can discretise the Neumann boundary condition
$\frac{\partial u}{\partial x_2} = 0$ at the nodes
$\boldsymbol{x}_{i,0} \;$, $i=0,\ldots,m$, by using the finite
difference approximation $$\frac{u_{i,0} - u_{i,-1}}{h} \ = \ 0 \ .$$
Therefore the second derivative at $\boldsymbol{x}_{i,0}$ becomes
$$- \ \frac{- u_{i,0} + u_{i,1}}{h^2} \qquad \text{instead of} \qquad
- \ \frac{u_{i,-1} - 2 u_{i,0} + u_{i,1}}{h^2}.$$ We can proceed in a
similar way at the other Neumann boundary nodes
$\{\boldsymbol{x}_{i,m} : i=0,\ldots,m_\delta \}$,
$\{\boldsymbol{x}_{0,j} : j=0,\ldots,m \}$ and
$\{\boldsymbol{x}_{m,j} : j=0,\ldots,m-1 \}$,  e.g.
$$- \ \frac{u_{m-1,j} - 2 u_{m,j} + u_{m+1,j}}{h^2}  \qquad
\text{is replaced by} \qquad - \ \frac{u_{m-1,j} - u_{m,j}}{h^2} \ .$$
On the other hand, zero coefficients $u_{i,m}$ for
$i = m_\delta+1,\ldots,m$, corresponding to the Dirichlet boundary
nodes, are not stored at all.

We order the indices $(i,j)$ lexicographically again, i.e. the unknowns
$u_{i,j}$ are stored in a vector $\boldsymbol{U}$ in the order
$(0,0),(1,0),
\ldots,(m,0), (0,1),\ldots,$ $(m,1),\ldots,(0,m),\ldots,(m_\delta,m)$.
Then the above equations form a system of
$$n \ := \ m\cdot(m+1) \ + \ m_\delta \ + \ 1$$ nonlinear equations in
$n$ unknowns (for simplicity we will not write down the dependency of
$\mathbf{G}$ on $\lambda$ and $\beta$ from now on): $$\label{model}
\mathbf{F}(\mathbf{U}) \ := \ A_\delta \, \mathbf{U}\ - \
\mathbf{G}(\mathbf{U}) \ = \ 0\vspace{1ex}$$ where
$\boldsymbol{U} := (u_{0,0},u_{1,0},\ldots,u_{m_\delta,m})^T \in \mathbb{R}^n$,
$\mathbf{G}(\mathbf{U}) := (g(u_{0,0}),\ldots,g(u_{m_\delta,m}))^T \in \mathbb{R}^n$,
and $A_\delta$ is the $n \times n$ matrix
$$A_\delta \ := \ \frac{1}{h^2} \left(
\begin{array}{cccccc}
B & -I \\
-I & C & -I \\
& -I & C & \ddots \\
& & \ddots & \ddots & \ddots \\
& & & \ddots & C & -I_\delta^T \\[0.5ex]
& & & & -I_\delta & B_\delta
\end{array}
\right)$$ with $(m+1)\times(m+1)$ matrices $$\arraycolsep4pt
B \ := \ \left(
\begin{array}{cccccc}
2 & -1 \\
-1 & 3 & -1 \\
& -1 & 3 & \ddots \\
& & \ddots & \ddots & \ddots \\
& & & \ddots & 3 & -1 \\[0.5ex]
& & & & -1 & 2
\end{array}
\right), \quad C \ := \ \left(
\begin{array}{cccccc}
3 & -1 \\
-1 & 4 & -1 \\
& -1 & 4 & \ddots \\
& & \ddots & \ddots & \ddots \\
& & & \ddots & 4 & -1 \\[0.5ex]
& & & & -1 & 3
\end{array}
\right) \ .$$ The $(m_\delta+1) \times m$ matrix $I_\delta$ is given by
$[I \ 0]$, where $I$ is the identity matrix, and the
$(m_\delta+1) \times (m_\delta+1)$ matrix $B_\delta$ consists of the
first $m_\delta+1$ rows and columns of $B$.

Now the Jacobian matrix $\mathbf{F}'(\mathbf{U})$ can be written
similarly to that in the first assignment, $$\label{eq:jac}
 \mathbf{F}'(\mathbf{U}) = A_\delta - \mathbf{G}'(\mathbf{U}),$$ where
$\mathbf{G}'(\mathbf{U})$ is a diagonal matrix with the values
$g'(u_{i,j})$ on the diagonal.

Inexact Newton Method {#inexact-newton-method .unnumbered}
---------------------

To solve the system of nonlinear equations
([\[model\]](#model){reference-type="ref" reference="model"}) we will
use the *Inexact Newton Method*:

Choose an *initial guess* $\boldsymbol{U}_0 \in \mathbb{R}^n$, a
*nonlinear tolerance* $\tau >0$ and *linear tolerances*
$\varepsilon_k>0$. Compute
$\mathbf{r}_0 = \mathbf{F}(\boldsymbol{U}_0)$.
($\|\mathbf{r}_k \|\leq \tau$) **exit** Approximately solve
$\mathbf{F}'(\boldsymbol{U}_k) \mathbf{s}_k = - \mathbf{r}_k$ for the
Newton step $\mathbf{s}_k$ up to
$\left\|\boldsymbol{r}^{\text{lin}}_{\ell}\right\| < \varepsilon_k \left\|\boldsymbol{r}^{\text{lin}}_0\right\|$,
as defined in [\[stop\]](#stop){reference-type="eqref"
reference="stop"}. Solution update:
$\mathbf{U}_{k+1} = \mathbf{U}_k + \mathbf{s}_k$ Residual update:
$\mathbf{r}_{k+1} = \mathbf{F}(\mathbf{U}_{k+1})$.

To solve the linear systems in Line 5 for each $k$, we will use the
Conjugate Gradient (CG) method with initial guess $\boldsymbol{0}$ and
stopping criterion $$\label{stop}
\|\boldsymbol{r}^{\text{lin}}_\ell \|\ \leq \
\varepsilon_k \|\boldsymbol{r}^{\text{lin}}_0 \|,$$ where
$\boldsymbol{r}^{\text{lin}}_\ell$ denotes the [linear]{.underline}
residual after $\ell$ CG iterations. For an arbitrary linear system
$$\label{eq:lin}
 A \boldsymbol{x} = \boldsymbol{b},$$ the linear residual is defined as
$\boldsymbol{r}^{\text{lin}}_{\ell}:=\boldsymbol{b} - A\boldsymbol{x}_{\ell}$,
where $\boldsymbol{x}_{\ell}$ is an approximate solution after $\ell$
iterations. The system [\[eq:lin\]](#eq:lin){reference-type="eqref"
reference="eq:lin"} with a symmetric positive definite matrix $A$ can be
solved by the **Conjugate Gradient Algorithm**:

Set $\boldsymbol{x}_0=\mathbf{0}$,
$\boldsymbol{r}^{\text{lin}}_{0}=\boldsymbol{b}$,
$\boldsymbol{q}_0=\boldsymbol{b}$.
$\alpha_\ell = ((\boldsymbol{r}_\ell^{\text{lin}})^T\boldsymbol{r}_\ell^{\text{lin}})/(\boldsymbol{q}_\ell^T A \boldsymbol{q}_\ell) \qquad\qquad\qquad\qquad$
$\boldsymbol{x}_{\ell+1} = \boldsymbol{x}_\ell + \alpha_\ell \boldsymbol{q}_\ell$
$\boldsymbol{r}_{\ell+1}^{\text{lin}} = \boldsymbol{r}_\ell^{\text{lin}} - \alpha_\ell A \boldsymbol{q}_\ell$
($\| \boldsymbol{r}^{\text{lin}}_{\ell+1} \| < \varepsilon_k \| \boldsymbol{r}^{\text{lin}}_0 \|$)
**exit**
$\beta_{\ell+1} = ((\boldsymbol{r}^{\text{lin}}_{\ell+1})^T\boldsymbol{r}^{\text{lin}}_{\ell+1})/((\boldsymbol{r}^{\text{lin}}_\ell)^T\boldsymbol{r}^{\text{lin}}_\ell)$
$\boldsymbol{q}_{\ell+1} = \boldsymbol{r}^{\text{lin}}_\ell + \beta_{\ell+1}
\boldsymbol{q}_\ell$

Note that the dot product of any two vectors
$\boldsymbol{a},\boldsymbol{b}\in\mathbb{R}^n$ is defined as
$\boldsymbol{a}^T \boldsymbol{b}:=\sum_{j=1}^n a_jb_j$ and hence the
quotients $\alpha_\ell$ and $\beta_{\ell+1}$ are real numbers.

Two reasonable choices for the linear tolerances $\varepsilon_k$ in the
stopping criterion ([\[stop\]](#stop){reference-type="ref"
reference="stop"}) are:

1.  fixed tolerance, e.g. $\varepsilon_k = 10^{-5}$

2.  variable tolerance: $\varepsilon_k = \min\left( \bar{\varepsilon} ,
    \gamma\cdot\dfrac{\|\mathbf{F}(\mathbf{U}_{k})\|^2}{\|\mathbf{F}(\mathbf{U}_{k-1})\|^2}\right)$
    with  $\varepsilon_0 = \bar{\varepsilon} = 0.1$  and  $\gamma =
    0.9$.

In this assignment only the second option is investigated. For more
details on the inexact Newton method and on the iterative solution of
systems of nonlinear equations in general see Kelley [@Ke:95].

The Code {#the-code .unnumbered}
========

main.f90 -

:   main program

header.f90 -

:   module containing the (parallel) sparse data structures (as in
    lectures)

laplace.f90 -

:   contains the routine Laplace() which sets up the matrix $A_{0.5}$ in
    compressed row format

save\_solution.f90 -

:   contains a routine that writes the solution vector $\boldsymbol{U}$
    to a file for postprocessing in Matlab (works only in sequential
    mode, i.e. for one processor)

maximum.f90 -

:   contains a function Maximum() which calculates the maximum of a
    distributed vector

newton.f90 -

:   a skeleton for the Newton() routine

func.f90 -

:   a skeleton for the Func() routine

jacobian.f90 -

:   a skeleton for the Jacobian() routine

cg.f90 -

:   contains the routine CG() which implements the CG method for the
    iterative solution of sparse linear systems in compressed row
    storage format. For a description of the arguments that CG() takes
    and of their intent see the comments at the beginning of cg.f90

matmult.f90 -

:   contains the routine Mat\_Mult() which implements the *sequential*
    sparse matrix vector product

vecdot.f90 -

:   a skeleton for the Vec\_Dot() function.

Makefile -

:   Makefile to compile and link all the necessary files on Balena

visualise.m -

:   Matlab function to visualise $u$ as a 3D surface plot

10 Adler J., Thermal--explosion theory for a slab with partial
insulation, *Combustion and Flame* **50**, 1983, pp. 1--7.  (Library:
PER66)

Greenway P. and Spence A., Numerical calculation of critical points for
a slab with partial insulation, *Combustion and Flame* **62**, 1985, pp.
141--156.  (Library: PER66)

Kelley CT., *Iterative Methods for Linear and Nonlinear Equations*,
SIAM, Philadelphia, 1995.  (Library: 512.978KEL)
