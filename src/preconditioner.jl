abstract type Preconditioner end
apply(p::Preconditioner, r) = throw("notimplemented")

@doc raw"""
    struct TracePreconditioner <: Preconditioner

Simple preconditioner for the interface problem. Uses the sum of the trace 
``tr(\cdot)`` of local stiffness matrices ``{K}_{bb}^{s}`` on the interface 
``I`` with interface dofs denoted by ``b`` and subdomains by ``s`` for ``N_s`` 
subdomains.

``P_{I}^{-1} = \sum_{s=1}^{N_s} tr(K_{bb}^{s}).``

# References
[1] : Equation (24) from Farhat, C. A Method of Finite Element Tearing and 
Interconnecting and its Parallel Solution Algorithm. International Journal for 
Numerical Methods in Engineering. vol 32, 1205-1227. 1991. 
"""
struct TracePreconditioner <: Preconditioner end

@doc raw"""
    struct DirichletPreconditioner <: Preconditioner

Dirichlet preconditioner that uses a Schur-complement to precondition the 
interface problem. This is the most accurate but also most expensive 
preconditioner.

The exact definition for the preconditioner ``D_{I}^{-1}`` for ``N_s`` 
subdomains is given below for interior dofs ``i``, interface dofs ``b``, 
subdomain ``s``, interface ``I``, and connectivity operator ``B``

``D_{I}^{-1} = \sum_{s=1}^{N_s} B^{s} \left[\begin{array}{c} 0 & 0 \\ 0 & S_{bb}^{s} \end{array}\right] B^{s^{T}},``

where the Schur complement ``S_{bb}^{s} = K_{bb}^{s} - K_{ib}^{s^T}K_{ii}^{s^{-1}}K_{ib}^{s}``.

# References
[1] : Equation (38) and (43) from Farhat, C., Mandel, J., and Roux, F.X. Optimal 
convergence properties of the FETI domain decomposition method. Computer Methods 
in Applied Mechanics and Engineering. vol 115, 365-385. 1994.
"""
struct DirichletPreconditioner <: Preconditioner end  

@doc raw"""
    struct LumpedPreconditioner <: Preconditioner

Preconditions the dual-interface problem with an economical preconditioner
only involving the the stiffness matrix ``K_{j}^{I}`` for a subdomain ``j``
at the interface ``I``.

The exact definition for the preconditioner ``L_{I}^{-1}``for ``N_s`` subdomains 
is given below for interface dofs ``b``, subdomain ``s``, and connectivity 
operator ``B``

``L_{I}^{-1} = \sum_{s=1}^{N_s} B^{s} \left[\begin{array}{c} 0 & 0 \\ 0 & K_{bb}^{s} \end{array}\right] B^{s^{T}}.``

# References
[1] : Equation (44) from Farhat, C., Mandel, J., and Roux, F.X. Optimal 
convergence properties of the FETI domain decomposition method. Computer Methods 
in Applied Mechanics and Engineering. vol 115, 365-385. 1994.
"""
struct LumpedPreconditioner <: Preconditioner end 
