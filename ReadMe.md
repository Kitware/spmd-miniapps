# SPMD Mini Applications

This is a collection of mini visualization applications, implemented
to analyze the suitability of vectorization for some popular
visualization algorithms.
The goal of this project is to evaluate the benefits and ease of implementation
of vectorized algorithm and to evaluate different tools that enable
vectorization.
In this project, we are primarily targeting SIMD capabilities of Intel CPUs,
especially the Xeon Phi.
Currently, we have implementations in
[ISPC](https://ispc.github.io/index.html).
The results of this experiment will help is in improving the performance of the
algorithms in VTK/ParaView.

# Algorithms

We have implemented the following algorithms:
* Structured Grid Contouring (Marching Cubes)
* Unstructured Gird Contouring (Marching Tetrahedra)
* 2D Poisson equation solver using SOR (successive over-relaxation)

Running each mini-app without any arguments will display its usage.
Each algorithm has a basic straight-forward scalar implementation and some
scalar and vectorized implementations, which try some variation in data
structures and data layouts.

## 2D Poisson Solver
This is an easy to vectorize algorithm. It was implemented to compare it
against the more complex algorithms (Marching Cubes and Marching Tetrahedra).
The application takes a grayscale pgm image as input, computes the divergence
of its gradients and runs SOR based Poisson solver on the divergence to get
back the original image.
Currently it does not produce the exact original image as output since we have
not implemented the initial values for the SOR and it may take many iterations
for the solver to converge.
Still, it serves its purpose as an easily vectorized implementation.

# Data Formats
All the applications work with legacy vtk file formats since readers for
these are quite simple to implement.
Datasets in any format can be easily converted to the legacy vtk format using
ParaView.
We have tested the applications with some of the datasets found here:
http://placid.nlm.nih.gov/community/21.
For tetrahedra mesh datasets, open any volumetric dataset in ParaView,
apply the Tetrahedralize filter and save the resulting dataset in the legacy
vtk format.

