abstract type Preconditioner end

"""
    struct DirichletPreconditioner <: Preconditioner

Dirichlet preconditioner that uses a Schur-complement to estimate interface 
forces. This is the most accurate but also most expensive preconditioner.
"""
struct DirichletPreconditioner <: Preconditioner end  

@doc raw"""
    struct TracePreconditioner <: Preconditioner

Preconditioner ``P_{I}^{-1} = \sum_{j=1}^{N_s} \mathbf{K}_{j}^{I}`` for 
computation of interface ``I`` constraints by the preconditioned conjugate 
projected gradient method for ``N_s`` subdomains with interface dofs of local 
stiffness matrices ``\mathbf{K}_{j}^{I}``.

# References
[1] : Equation (24) from Farhat, C. A Method of Finite Element Tearing and 
Interconnecting and its Parallel Solution Algorithm. vol 32, 1205-1227. 1991.  
"""
struct TracePreconditioner <: Preconditioner end
