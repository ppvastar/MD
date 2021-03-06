# Read Me
This directory includes the source code for the numerical calculations in manuscript http://arxiv.org/abs/1601.01402.
The file "Makefile" is used to compile the main source file "main.f90" on the Linux platform.

Note that two preprocessor directives appear in the Makefile:
    
MACRO+=-D__solve_topological_charge

MACRO+=-D__solve_eigenstates

The first directive is used to calculate the topological charge of the QD chain with varying Zeeman
energy and chemical potential [e.g., as shown by Fig. 2(a) in the manuscript]. The second directive is used to solve 
the eigenstates of the QD chain [e.g., as shown by Fig. 2(b) in the manuscript]. The directive does not take effect if 
it is commented by "#" at the left beginning.

The MKL Lapack library is required. 
