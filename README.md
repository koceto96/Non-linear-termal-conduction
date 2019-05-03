\\documentclass\[11pt,a4paper\]{article} \\usepackage{amsmath}
\\usepackage{amssymb} \\usepackage{epsfig} \\usepackage{comment}
\\usepackage{color} \\usepackage{algorithm,algpseudocode}
\\usepackage{stmaryrd}

\\newcommand{\\degree}{{}\^{\\circ}}
\\newcommand{\\DDt}\[1\]{\\frac{D\#1}{Dt}}
\\newcommand{\\ddz}\[1\]{\\frac{\\partial \#1}{\\partial z}}
\\renewcommand{\\vec}\[1\]{\\mathbf{\#1}} \\newcommand{\\dt}{\\Delta t}
\\newcommand{\\adt}{\\alpha\\dt} \\newcommand{\\omadt}{(1-\\alpha)\\dt}
\\newcommand{\\const}{\\mathrm{const}}
\\newcommand{\\refr}{{\\mathrm{ref}}}
\\newcommand{\\order}\[1\]{\\mathcal{O}(\#1)}
\\newcommand{\\Ukelvin}{\\operatorname{K}}
\\newcommand{\\Um}{\\operatorname{m}}
\\newcommand{\\Us}{\\operatorname{s}}
\\newcommand{\\Ukg}{\\operatorname{kg}}
\\newcommand{\\Rearth}{R\_{\\operatorname{earth}}}
\\newcommand{\\bg}{{\*0}} \\newcommand{\\Id}{\\operatorname{Id}}
\\def\\bF{\\mathbf{F}} \\def\\bU{\\mathbf{U}} \\def\\br{\\mathbf{r}}
\\def\\bs{\\mathbf{s}} \\def\\Vert{\|}

\\newcommand{\\cchange}\[1\]{{\\color{black}\#1}}
\\newcommand{\\cred}\[1\]{{\\color{red}\#1}}
\\newcommand{\\stilltodo}\[1\]{\\cred{\#1}}
\\newcommand{\\dune}{\\textsc{Dune}}
\\newcommand{\\alugrid}{\\textsc{ALUGrid}}
\\newcommand{\\uggrid}{\\textsc{UGGrid}}
\\newcommand{\\yaspgrid}{\\texttt{YaspGrid}}
\\newcommand{\\aquila}{\\textsc{Aquila}}
\\newcommand{\\openmpi}{\\textsc{OpenMPI}}
\\newcommand{\\ENDGame}{\\textsc{ENDGame}}
\\newcommand{\\newdynamics}{\\textsc{New Dynamics}}
\\newcommand{\\metoffice}{\\textsc{Met Office}}

\\pagestyle{plain} \\pagenumbering{arabic} \\raggedbottom \\topmargin1cm
\\headheight0pt \\headsep30pt \\topskip0pt \\textheight21.5cm
\\oddsidemargin0.6cm \\evensidemargin0.6cm \\textwidth15.5cm

\\setlength{\\oddsidemargin}{0.0cm}
\\setlength{\\evensidemargin}{-0.04cm} \\setlength{\\textwidth}{15.5cm}
\\setlength{\\topmargin}{-1.54cm} \\setlength{\\textheight}{24.9cm}
%22.2 \\renewcommand{\\hoffset}{0cm} \\renewcommand{\\voffset}{0cm}
\\newcommand{\\laplacehoriz}{\\Delta\^{(2\\operatorname{d})}}

\\usepackage{tikz}
\\usetikzlibrary{positioning,shapes,calc,backgrounds,patterns}
\\usepackage{pgfplots}

\\begin{document} %%%%%%%%%%%%%%%%%%%% \\section\*{Nonlinear Thermal
Conduction} %%%%%%%%%%%%%%%%%%%%%%%%% This project is about practical
aspects of solving sparse systems of nonlinear equations
\\begin{equation} \\label{nonlin} \\bF(\\boldsymbol{U})  =
 \\boldsymbol{0} \\end{equation} using the inexact Newton method. In
particular we will use parallel Newton\--CG and look at a case study on
nonlinear thermal combustion in a self\--heating medium in a partially
insulated square domain. This is the same model as the one considered in
the first project. This problem arises when investigating critical
parameter values for underground repositories of self\--heating waste
which are partially covered by buildings. Beyond certain critical
parameter values the steady state solutions may become very large or
even unbounded leading to an explosion in the medium. For more
information see Greenway & Spence \\cite{GrSp:85} and Adler
\\cite{Ad:83}.

\\subsection\*{Model formulation} We are concerned with finding
solutions of the (dimensionless) nonlinear thermal conduction equation
\\begin{equation} \\label{pde} \\mathcal{F}(u)\[\\lambda,\\beta\]  :=  -
 \\frac{\\partial\^2 u}{\\partial x\_1\^2}  -
\\frac{\\partial\^2 u}{\\partial x\_2\^2}  -  g(u)\[\\lambda,\\beta\]  =
 0 \\qquad \\text{in } \\Omega := \[0,1\]\^2, \\end{equation} subject to
the following mixed boundary conditions on \$\\partial \\Omega\$:

\\begin{minipage}\[b\]{0.40\\linewidth} \\begin{equation*} % \\label{bc}
\\hspace*{-2cm} \\begin{array}{c} \\text{Dirichlet}
\\left{\\begin{array}{rcll} u & = & 0 & \\quad \\text{if }  x\_1 \>
\\delta  \\text{ and }  x\_2 = 1 \\end{array}\\right.\\\[1ex\]
\\text{Neumann} \\left{\\begin{array}{rcll} \\displaystyle
\\frac{\\partial u}{\\partial x\_1} & = & 0 & \\quad \\text{if }  x\_1 =
0  \\text{ or }  x\_1 = 1\\\[2ex\] \\displaystyle \\frac{\\partial
u}{\\partial x\_2} & = & 0 & \\quad \\text{otherwise} \\end{array}
\\right. \\end{array} \\end{equation*} \\end{minipage} \\hfill
\\begin{minipage}\[c\]{0.46\\linewidth} \\begin{tikzpicture}
\\node\[rectangle,draw=black,minimum width=3.2cm,minimum
height=3.2cm,inner sep=0cm\] (omega) at (0,0)
{\$\\mathcal{F}(u)\[\\lambda,\\beta\]=0\$}; \\node\[anchor=south east\]
at (omega.south west) {\$0\$}; \\node\[anchor=north west\] at
(omega.south west) {\$0\$}; \\node\[anchor=north\] at
(\$0.3*(omega.south west)+0.7*(omega.south east)\$) {\$\\delta\$};
\\draw\[-,black\] (\$0.3*(omega.south west)+0.7*(omega.south
east)+(0,-0.1)\$) \-- (\$0.3*(omega.south west)+0.7*(omega.south
east)+(0,0.1)\$); \\node\[anchor=center,ellipse,draw=black,minimum
width=5pt,minimum height=5pt,inner sep=0pt\] at (\$0.3*(omega.north
west)+0.7*(omega.north east)\$) {}; \\node\[anchor=north east\] at
(omega.south east) {\$1\$}; \\node\[anchor=north east\] at (omega.north
west) {\$1\$}; \\node\[anchor=north,inner sep=16pt\] at (omega.south)
{\$\\frac{\\partial u}{\\partial x\_2}=0\$}; \\node\[anchor=east,inner
sep=8pt\] at (omega.west) {\$\\frac{\\partial u}{\\partial x\_1}=0\$};
\\node\[anchor=west,inner sep=8pt\] at (omega.east) {\$\\frac{\\partial
u}{\\partial x\_1}=0\$}; \\node\[anchor=south\] at (\$0.7*(omega.north
west)+0.3*(omega.north east)\$) {\$\\frac{\\partial u}{\\partial
x\_2}=0\$}; \\node\[anchor=south,inner sep=9pt\] at (\$0.1*(omega.north
west)+0.9\*(omega.north east)\$) {\$u=0\$}; \\draw\[-\>,black\]
(-1.6cm,-1.6cm) \-- (-1.6cm,2.3cm) node \[at end,anchor=north east\]
{\$x\_2\$}; \\draw\[-\>,black\] (-1.6cm,-1.6cm) \-- (2.5cm,-1.6cm) node
\[at end,anchor=north east\] {\$x\_1\$}; \\end{tikzpicture}
\\end{minipage}

\\noindent for \$0 \\le \\delta \< 1\$. Note that these (more realistic)
boundary conditions differ from those from Assignment\~1. The nonlinear
functional \\begin{equation} \\label{arrhen} g(u)\[\\lambda,\\beta\]  :=
 \\lambda ; \\exp\\left( \\frac{u}{1+\\beta u} \\right) \\end{equation}
is the Arrhenius reaction rate which depends on the two parameters
\$\\lambda\$ and \$\\beta\$.

The relationship between the dimensionless variables in (\\ref{pde}),
(\\ref{arrhen}) and the physical quantities is given by Adler
\\cite{Ad:83}. We note that \$u\$ is a dimensionless temperature excess
and \$\\lambda\$ the Frank\--Kamenetskii parameter; \$\\beta\$ is the
dimensionless activation energy and \$\\delta\$ the dimensionless
half\--width of an insulating strip.

\\subsection\*{Finite Difference Discretisation}

As for the linear Poisson equation, suppose that \$\\Omega\$ is
subdivided into small equal squares of size \$h \\times h\$ with \$h :=
1/m\$. The interior and boundary nodes of the mesh are given by
\$\\boldsymbol{x}*{i,j} := (ih,jh)\$, with \$i,j = 0,\\ldots,m\$. Let
\$m*\\delta := \[\\delta m\] \\in \\mathbb{N}\$ be the integer part of
\$\\delta m\$ (i.e. largest integer \$m\_{\\delta}\\le \\delta m\$).
Then \$m\_\\delta h \\le \\delta\$ and \$(m\_\\delta + 1)h \> \\delta\$,
and so all the points \$\\boldsymbol{x}*{i,m}\$ with \$i \> m*\\delta\$
lie on the Dirichlet boundary.

The finite difference approach to the above problem is now to seek
approximations \$u\_{i,j}\$ to \$u(\\boldsymbol{x}*{i,j})\$ such that \[
- \\frac{u*{i-1,j} - 2 u\_{i,j} + u\_{i+1,j}}{h\^2}  -
\\frac{u\_{i,j-1} - 2 u\_{i,j} + u\_{i,j+1}}{h\^2}  -
g(u\_{i,j})  =  0. \] We can discretise the Neumann boundary condition
\$\\frac{\\partial u}{\\partial x\_2} = 0\$ at the nodes\\linebreak
\$\\boldsymbol{x}*{i,0} ;\$, \$i=0,\\ldots,m\$, by using the finite
difference approximation \[ \\frac{u*{i,0} - u\_{i,-1}}{h}  =  0  . \]
Therefore the second derivative at \$\\boldsymbol{x}*{i,0}\$ becomes \[
-  \\frac{- u*{i,0} + u\_{i,1}}{h\^2} \\qquad \\text{instead of} \\qquad
-  \\frac{u\_{i,-1} - 2 u\_{i,0} + u\_{i,1}}{h\^2}. \] We can proceed in
a similar way at the other Neumann boundary nodes
\${\\boldsymbol{x}*{i,m} : i=0,\\ldots,m*\\delta }\$,
\${\\boldsymbol{x}*{0,j} : j=0,\\ldots,m }\$ and
\${\\boldsymbol{x}*{m,j} : j=0,\\ldots,m-1 }\$,  e.g. \[ -
 \\frac{u\_{m-1,j} - 2 u\_{m,j} + u\_{m+1,j}}{h\^2} \\qquad \\text{is
replaced by} \\qquad -  \\frac{u\_{m-1,j} - u\_{m,j}}{h\^2}  . \] On the
other hand, zero coefficients \$u\_{i,m}\$ for \$i =
m\_\\delta+1,\\ldots,m\$, corresponding to the Dirichlet boundary nodes,
are not stored at all.

We order the indices \$(i,j)\$ lexicographically again, i.e. the
unknowns \$u\_{i,j}\$ are stored in a vector \$\\boldsymbol{U}\$ in the
order \$(0,0),(1,0), \\ldots,(m,0), (0,1),\\ldots,\$
\$(m,1),\\ldots,(0,m),\\ldots,(m\_\\delta,m)\$. Then the above equations
form a system of \[ n  :=  m\\cdot(m+1)  +  m\_\\delta  +  1 \]
nonlinear equations in \$n\$ unknowns (for simplicity we will not write
down the dependency of \$\\mathbf{G}\$ on \$\\lambda\$ and \$\\beta\$
from now on):\\vspace{1ex} \\begin{equation} \\label{model} \\bF(\\bU)
 :=  A\_\\delta , \\bU  -
\\mathbf{G}(\\bU)  =  0\\vspace{1ex} \\end{equation} where
\$\\boldsymbol{U} := (u\_{0,0},u\_{1,0},\\ldots,u\_{m\_\\delta,m})\^T
\\in \\mathbb{R}\^n\$, \$\\mathbf{G}(\\bU) :=
(g(u\_{0,0}),\\ldots,g(u\_{m\_\\delta,m}))\^T \\in \\mathbb{R}\^n\$, and
\$A\_\\delta\$ is the \$n \\times n\$ matrix \[ A\_\\delta  :=
 \\frac{1}{h\^2} \\left( \\begin{array}{cccccc} B & -I \\ -I & C & -I \\
& -I & C & \\ddots \\ & & \\ddots & \\ddots & \\ddots \\ & & & \\ddots &
C & -I\_\\delta\^T \\\[0.5ex\] & & & & -I\_\\delta & B\_\\delta
\\end{array} \\right) \] with \$(m+1)\\times(m+1)\$ matrices \[
\\arraycolsep4pt B  :=  \\left( \\begin{array}{cccccc} 2 & -1 \\ -1 & 3
& -1 \\ & -1 & 3 & \\ddots \\ & & \\ddots & \\ddots & \\ddots \\ & & &
\\ddots & 3 & -1 \\\[0.5ex\] & & & & -1 & 2 \\end{array} \\right),
\\quad C  :=  \\left( \\begin{array}{cccccc} 3 & -1 \\ -1 & 4 & -1 \\ &
-1 & 4 & \\ddots \\ & & \\ddots & \\ddots & \\ddots \\ & & & \\ddots & 4
& -1 \\\[0.5ex\] & & & & -1 & 3 \\end{array} \\right)  . \] The
\$(m\_\\delta+1) \\times m\$ matrix \$I\_\\delta\$ is given by \$\[I
 0\]\$, where \$I\$ is the identity matrix, and the \$(m\_\\delta+1)
\\times (m\_\\delta+1)\$ matrix \$B\_\\delta\$ consists of the first
\$m\_\\delta+1\$ rows and columns of \$B\$.

Now the Jacobian matrix \$\\bF\'(\\bU)\$ can be written similarly to
that in the first assignment, \\begin{equation}\\label{eq:jac}
\\bF\'(\\bU) = A\_\\delta - \\mathbf{G}\'(\\bU), \\end{equation} where
\$\\mathbf{G}\'(\\bU)\$ is a diagonal matrix with the values
\$g\'(u\_{i,j})\$ on the diagonal.

\\newpage \\subsection\*{Inexact Newton Method}

To solve the system of nonlinear equations (\\ref{model}) we will use
the {\\em Inexact Newton Method}: \\begin{algorithm}\[h!\]
\\caption{Inexact Newton} \\label{alg:new} \\begin{algorithmic}\[1\]
\\State Choose an {\\em initial guess} \$\\boldsymbol{U}\_0 \\in
\\mathbb{R}\^n\$, a {\\em nonlinear tolerance} \$\\tau \>0\$ and {\\em
linear tolerances} \$\\varepsilon\_k\>0\$. \\State Compute \$\\br\_0 =
\\bF(\\boldsymbol{U}*0)\$. \\For{\$k=0,1,2,\\dots,k*{\\max}-1\$} \\State
{\\bf If} (\$\\Vert \\br\_k \\Vert \\leq \\tau\$) {\\bf exit} \\State
Approximately solve \$\\bF\'(\\boldsymbol{U}\_k) \\bs\_k = - \\br\_k\$
for the Newton step \$\\bs\_k\$ \\State \\hspace\*{0.17\\linewidth} up
to \$\\left\|\\boldsymbol{r}\^{\\text{lin}}*{\\ell}\\right\| \<
\\varepsilon\_k \\left\|\\boldsymbol{r}\^{\\text{lin}}*0\\right\|\$, as
defined in \\eqref{stop}. \\State Solution update: \$\\bU*{k+1} =
\\bU\_k + \\bs\_k\$ \\State Residual update: \$\\br*{k+1} =
\\bF(\\bU\_{k+1})\$. \\EndFor \\end{algorithmic} \\end{algorithm}

To solve the linear systems in Line 5 for each \$k\$, we will use the
Conjugate Gradient (CG) method with initial guess \$\\boldsymbol{0}\$
and stopping criterion \\begin{equation} \\label{stop} \\Vert
\\boldsymbol{r}\^{\\text{lin}}*\\ell \\Vert  \\leq
\\varepsilon\_k \\Vert \\boldsymbol{r}\^{\\text{lin}}*0 \\Vert,
\\end{equation} where \$\\boldsymbol{r}\^{\\text{lin}}*\\ell\$ denotes
the \\underline{linear} residual after \$\\ell\$ CG iterations. For an
arbitrary linear system \\begin{equation}\\label{eq:lin} A
\\boldsymbol{x} = \\boldsymbol{b}, \\end{equation} the linear residual
is defined as \$\\boldsymbol{r}\^{\\text{lin}}*{\\ell}:=\\boldsymbol{b}
- A\\boldsymbol{x}*{\\ell}\$, where \$\\boldsymbol{x}*{\\ell}\$ is an
approximate solution after \$\\ell\$ iterations. The system
\\eqref{eq:lin} with a symmetric positive definite matrix \$A\$ can be
solved by the \\textbf{Conjugate Gradient Algorithm}:
\\begin{algorithm}\[h!\] \\caption{Conjugate gradient (CG) method}
\\label{alg:cg} \\begin{algorithmic}\[1\] \\State Set
\$\\boldsymbol{x}*0=\\mathbf{0}\$,
\$\\boldsymbol{r}\^{\\text{lin}}*{0}=\\boldsymbol{b}\$,
\$\\boldsymbol{q}*0=\\boldsymbol{b}\$.
\\For{\$\\ell=0,1,2,\\ldots,\\ell*{\\max}-1\$} \\State \$\\alpha\_\\ell
=
((\\boldsymbol{r}*\\ell^{\\text{lin}})^T\\boldsymbol{r}*\\ell\^{\\text{lin}})/(\\boldsymbol{q}*\\ell\^T
A \\boldsymbol{q}*\\ell) \\qquad\\qquad\\qquad\\qquad\$ \\Comment{search
parameter} \\State \$\\boldsymbol{x}*{\\ell+1} = \\boldsymbol{x}*\\ell +
\\alpha\_\\ell \\boldsymbol{q}*\\ell\$ \\Comment{update solution}
\\State \$\\boldsymbol{r}*{\\ell+1}\^{\\text{lin}} =
\\boldsymbol{r}*\\ell\^{\\text{lin}} - \\alpha*\\ell A
\\boldsymbol{q}*\\ell\$ \\Comment{update residual} \\State {\\bf If}
(\$\|
\\boldsymbol{r}^{\\text{lin}}*{\\ell+1}\ \|\ \<\ \\varepsilon\_k\ \|\ \\boldsymbol{r}\^{\\text{lin}}*0\ \|\$)\ {\\bf\ exit}\ \\State\ \$\\beta*{\\ell+1}\ =\ ((\\boldsymbol{r}^{\\text{lin}}*{\\ell+1})^T\\boldsymbol{r}^{\\text{lin}}*{\\ell+1})/((\\boldsymbol{r}^{\\text{lin}}*\\ell)^T\\boldsymbol{r}\^{\\text{lin}}*\\ell)\$
\\State \$\\boldsymbol{q}*{\\ell+1} =
\\boldsymbol{r}\^{\\text{lin}}*\\ell + \\beta\_{\\ell+1}
\\boldsymbol{q}\_\\ell\$ \\Comment{update search direction} \\EndFor
\\end{algorithmic} \\end{algorithm}

\\noindent Note that the dot product of any two vectors
\$\\boldsymbol{a},\\boldsymbol{b}\\in\\mathbb{R}\^n\$ is defined as
\$\\boldsymbol{a}\^T \\boldsymbol{b}:=\\sum\_{j=1}\^n a\_jb\_j\$ and
hence the quotients \$\\alpha\_\\ell\$ and \$\\beta\_{\\ell+1}\$ are
real numbers.

Two reasonable choices for the linear tolerances \$\\varepsilon\_k\$ in
the stopping criterion (\\ref{stop}) are: \\begin{enumerate} \\item
fixed tolerance, e.g. \$\\varepsilon\_k = 10\^{-5}\$ \\item variable
tolerance: \$ \\varepsilon\_k = \\min\\left( \\bar{\\varepsilon} ,
\\gamma\\cdot\\dfrac{\\Vert\\bF(\\bU\_{k})\\Vert^2}{\\Vert\\bF(\\bU\_{k-1})\\Vert^2}\\right)
\$ with  \$\\varepsilon\_0 = \\bar{\\varepsilon} = 0.1\$  and  \$\\gamma
= 0.9\$. \\end{enumerate} In this assignment only the second option is
investigated. For more details on the inexact Newton method and on the
iterative solution of systems of nonlinear equations in general see
Kelley \\cite{Ke:95}.

\\newpage

\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \\section\*{The Code}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\\begin{description} \\item\[{\\tt main.f90} -\] main
program\\vspace{-1.5ex} \\item\[{\\tt header.f90} -\] module containing
the (parallel) sparse data structures (as in lectures)\\vspace{-1.5ex}
\\item\[{\\tt laplace.f90} -\] contains the routine {\\tt Laplace()}
which sets up the matrix \$A\_{0.5}\$ in compressed row
format\\vspace{-1.5ex} \\item\[{\\tt save\_solution.f90} -\] contains a
routine that writes the solution vector \$\\boldsymbol{U}\$ to a file
for postprocessing in {\\tt Matlab} (works only in sequential mode, i.e.
for one processor)\\vspace{-1.5ex} \\item\[{\\tt maximum.f90} -\]
contains a function {\\tt Maximum()} which calculates the maximum of a
distributed vector\\vspace{-1.5ex} \\item\[{\\tt newton.f90} -\] a
skeleton for the {\\tt Newton()} routine \\vspace{-1.5ex} \\item\[{\\tt
func.f90} -\] a skeleton for the {\\tt Func()} routine \\vspace{-1.5ex}
\\item\[{\\tt jacobian.f90} -\] a skeleton for the {\\tt Jacobian()}
routine \\vspace{-1.5ex} \\item\[{\\tt cg.f90} -\] contains the routine
{\\tt CG()} which implements the CG method for the iterative solution of
sparse linear systems in compressed row storage format. For a
description of the arguments that {\\tt CG()} takes and of their {\\tt
intent} see the comments at the beginning of {\\tt cg.f90}
\\vspace{-1.5ex} \\item\[{\\tt matmult.f90} -\] contains the routine
{\\tt Mat\_Mult()} which implements the \\emph{sequential} sparse matrix
vector product\\vspace{-1.5ex} \\item\[{\\tt vecdot.f90} -\] a skeleton
for the {\\tt Vec\_Dot()} function. \\vspace{-1.5ex} \\item\[{\\tt
Makefile} -\] Makefile to compile and link all the necessary files on
{\\tt Balena}\\vspace{-1.5ex} \\item\[{\\tt visualise.m} -\] Matlab
function to visualise \$u\$ as a 3D surface plot \\end{description}

\\begin{thebibliography}{10} \\bibitem{Ad:83} Adler J.,
Thermal\--explosion theory for a slab with partial insulation, {\\em
Combustion and Flame} {\\bf 50}, 1983, pp. 1\--7.  (Library: {\\tt
PER66})

\\bibitem{GrSp:85} Greenway P. and Spence A., Numerical calculation of
critical points for a slab with partial insulation, {\\em Combustion and
Flame} {\\bf 62}, 1985, pp. 141\--156.  (Library: {\\tt PER66})

\\bibitem{Ke:95} Kelley CT., {\\em Iterative Methods for Linear and
Nonlinear Equations}, SIAM, Philadelphia, 1995.  (Library: {\\tt
512.978KEL})

\\end{thebibliography}

\\end{document}
