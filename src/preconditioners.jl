abstract type Preconditioner end

"""
    struct DirichletPreconditioner <: Preconditioner

Dirichlet preconditioner that uses a Schur-complement to estimate interface 
forces. This is the most accurate but also most expensive preconditioner.
"""
struct DirichletPreconditioner <: Preconditioner

end  
