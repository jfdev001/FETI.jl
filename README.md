# FETI

<!---[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jfdev001.github.io/FETI.jl/stable/)--->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jfdev001.github.io/FETI.jl/dev/)
[![Build Status](https://github.com/jfdev001/FETI.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jfdev001/FETI.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

A prototype finite element tearing and interconnecting (FETI) solver. Serial algorithm will be implemented first, and possibly the parallel version will be included. A parallel version with a data oriented approach (e.g., using PartitionedArrays.jl) could be more simple than using MPI.

The pseudocode for the FETI solver itself:

```julia
# args could just be workspace...
# TODO: How is F_I not explicitly formed... perhaps see section 6 with parallel solve,
# check AMFeti, and check Krylov.jl or IterativeSolvers.jl for how the A*s operation in 
# in the pcg is done in a "matrix free" setting.
# TODO: How is r0 formed for the PCPG for N_s subdomains? 
function feti_solver(LocalProblems, Connector, InterfacePreconditioner, InterfaceSolver)
    lambda = initial_lambda(LocalProblems, Connector) # Appendix 2, equation (37)
    F_I = compute_F_I(LocalProblems, Connector) # eq(14), but sec. 5 not explicit assembled
    alpha = compute_alpha(lambda, F_I, LocalProblems, Connector) # pg. 1213
    projection_P = projection(LocalProblems) # pg. 1213
    set_projection_P(InterfaceSolver, projection_P)
    set_preconditioner(InterfaceSolver, InterfacePreconditioner)
    while !converged 
        # Direct solves on local problems in parallel
        # TODO: sec 6, 
        # "...all local finite element computations can be performed in parallel. 
        # [Including]... factoring Kj"... but how can Kj be factored if it is singular
        # on floating subdomains...
        for j, LocalProblem in enumerate(LocalProblems)
            # equation (8) and (10)
            Rj = nullspace(LocalProblem) # pg. 1213
            Kj_pinv = pinv(LocalProblem)
            fj = rhs(LocalProblem)
            BjT = sum(transpose.(connectivity(j, Connector)))
            uj = Kj_pinv*(fj + BjT) + Rj*alpha
        end
    
        # parallel'ish interface solve
        set_lambda(InterfaceSolver, lambda)
        lambda = solve(InterfaceSolver, LocalProblems, Connector) 
        alpha = compute_alpha(lambda, F_I, LocalProblems, Connector)
    end 
end 
```

Initial pseudocode below with slightly more details and less references to equations.

```julia
## "Verbose" pseudocode
# Setup local subdomain problems by using a partitioned mesh 
mesh = make_mesh()
partitions = partition_mesh(mesh) 
local_problems = []
for part_i, partition in enumerate(partitions): 
    ki, bi <- assemble_local_problem(partition)
    push!(local_problems, LocalProblem(ki, bi))
end 

# Setup the connectivity info between partitions so that it is known which
# dofs are on the interface and will be subject to continuity constraints 
# ... should essentially tell the indices of different interface vs interior 
# dofs   
connectivity_info = ConnectivityInfo()
for part_i in 1:length(partitions)
    for part_j in i+1:length(partitions) # all other partitions
        # get interface connectivity information 
        boundary_i = get_boundary(partitions[part_i])
        boundary_j = get_boundary(partitions[part_j])
        interface = intersection(boundary_i, boundary_j)
        update_connectivity_info!(connectivity_info, interface, part_i, part_j)
    end 
 
    # All dofs for part i not on the interface are interior dofs    
    all_dofs = get_dofs(partitions[part_i])
    interface_dofs = get_all_interface_dofs(connectivity_info, part_i)
    interior_dofs = all_dofs[all_dofs .∉ Ref(interface_dofs)]
end

# Given the delineation between interior and interface dofs, domain 
# decomposition methods are then concerned with handling the fact that
# local problems are not guaranteed to be non-singular (and thus admit 
# a unique solution)... so one performs manipulations to handle this...
...

# The continuity at the interface is in FETI enforced using langrange
# multipliers and solving a saddle point problem
... 
```

## Constraint Matrix

The Finite Element Tearing and Interconnecting (FETI) method is a domain decomposition technique used in solving large-scale systems of linear equations arising from finite element analysis. In FETI, the global problem is decomposed into smaller subproblems which are then solved in parallel, and a constraint matrix is used to enforce the continuity of the solution across subdomain interfaces.

### Constraint Matrix in FETI

The constraint matrix \( B \) in FETI is used to enforce the continuity conditions on the interfaces between subdomains. Each row of the matrix \( B \) represents a constraint that ensures that the degrees of freedom (DOFs) on the interface between two subdomains are equal.

### Example

Consider a simple example with two subdomains \( \Omega_1 \) and \( \Omega_2 \) that share a common interface. Let's denote:
- \( u_1 \) as the vector of DOFs in \( \Omega_1 \),
- \( u_2 \) as the vector of DOFs in \( \Omega_2 \).

Assume \( u_1 \) and \( u_2 \) have been partitioned such that the interface DOFs are the last entries in these vectors. Let \( u_{1,I} \) and \( u_{2,I} \) be the interface DOFs for subdomains \( \Omega_1 \) and \( \Omega_2 \), respectively.

The constraint that the interface DOFs must be equal can be written as:
\[ u_{1,I} = u_{2,I} \]

This can be represented in matrix form as:
\[ B \begin{pmatrix} u_1 \\ u_2 \end{pmatrix} = 0 \]

### Constructing the Constraint Matrix \( B \)

For simplicity, let's assume each subdomain has 3 DOFs, and the interface DOF is the 3rd DOF in each subdomain vector:
\[ u_1 = \begin{pmatrix} u_{1,1} \\ u_{1,2} \\ u_{1,I} \end{pmatrix}, \quad u_2 = \begin{pmatrix} u_{2,1} \\ u_{2,2} \\ u_{2,I} \end{pmatrix} \]

The constraint \( u_{1,I} = u_{2,I} \) can be written as:
\[ u_{1,3} - u_{2,3} = 0 \]

This can be represented using the constraint matrix \( B \):
\[ B = \begin{pmatrix} 0 & 0 & 1 & 0 & 0 & -1 \end{pmatrix} \]

So the constraint equation is:
\[ B \begin{pmatrix} u_{1,1} \\ u_{1,2} \\ u_{1,3} \\ u_{2,1} \\ u_{2,2} \\ u_{2,3} \end{pmatrix} = \begin{pmatrix} 0 & 0 & 1 & 0 & 0 & -1 \end{pmatrix} \begin{pmatrix} u_{1,1} \\ u_{1,2} \\ u_{1,3} \\ u_{2,1} \\ u_{2,2} \\ u_{2,3} \end{pmatrix} = 0 \]

### Generalization

In a more general case with multiple subdomains and multiple interface DOFs, the constraint matrix \( B \) will have multiple rows, each corresponding to a continuity constraint for an interface DOF pair.

For instance, if we have \( n \) subdomains, each with \( m \) DOFs, and a set of \( k \) interface constraints, the constraint matrix \( B \) will be a \( k \times (n \times m) \) matrix where each row enforces a particular interface DOF continuity constraint.

### Ensuring \( Bu = 0 \)

To ensure the system \( Bu = 0 \) is satisfied, we typically employ Lagrange multipliers. These multipliers are introduced to enforce the constraints during the minimization of the global energy function. The resulting system of equations, which includes both the original equations and the constraints, is solved simultaneously.

The augmented system of equations is:
\[ \begin{pmatrix} K & B^T \\ B & 0 \end{pmatrix} \begin{pmatrix} u \\ \lambda \end{pmatrix} = \begin{pmatrix} f \\ 0 \end{pmatrix} \]

where:
- \( K \) is the global stiffness matrix,
- \( B \) is the constraint matrix,
- \( u \) is the vector of unknown DOFs,
- \( \lambda \) is the vector of Lagrange multipliers,
- \( f \) is the global force vector.

By solving this augmented system, we ensure that the constraints \( Bu = 0 \) are satisfied, leading to a consistent and accurate solution across the subdomains.

# References

[1] : AMFeti. url: https://github.com/AppliedMechanics/AMfeti

[2] : Farhat, C. A Method of Finite Element Tearing and Interconnecting and its 
Parallel Solution Algorithm. International Journal for Numerical Methods in
Engineering. vol 32, 1205-1227. 1991. 

[3] : Farhat, C., Mandel, J., and Roux, F.X. Optimal convergence properties of 
the FETI domain decomposition method. Computer Methods in Applied Mechanics and
Engineering. vol 115, 365-385. 1994.

[4] : Lesoinne, M. and Person, K. An Efficient FETI Implementation on 
Distributed Shared Memory Machines with Independent Numbers of Subdomains and
Processors. American Mathematical Society. vol 218. 1998.

[5] : Dostal, Z., Horak, D., and Kucera, R. Total FETI - An Easier Implementable 
Variant of the FETI Method for Numerical Solution of Elliptic PDE. 
Communications in Numerical Methods in Engineering. vol 22 (12), 1155-1162. 2006
