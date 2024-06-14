# Combined connectivity, pcpg, and fem to solve a linear system of equations
# (assumes linear static problem) via FETI
export interface_solver, interface_preconditioner

"""
    struct ClassicalFETI

Classical finite element tearing and interconnecting (FETI) workspace.

NOTE: Could use lazy forwarding here to get accessors for composing classes

# References
[1] : Farhat, C. A Method of Finite Element Tearing and 
Interconnecting and its Parallel Solution Algorithm. International Journal for 
Numerical Methods in Engineering. vol 32, 1205-1227. 1991. 
"""
@kwdef struct ClassicalFETI{TLPS, TIS, TIP}
    local_problems::TLP
    interface_solver::TIS
    interface_preconditioner::TIP
end

local_problems(feti::ClassicalFETI) = feti.local_problems
interface_solver(feti::ClassicalFETI) = feti.interface_solver
interface_preconditioner(feti::ClassicalFETI) = feti.interface_preconditioner

"""
    ClassicalFETI(...)

TODO!

Constructs the `struct ClassicalFETI` by applying an `assemble` operation
on a partitioned (by METIS) `mesh` and constructing local problems to be solved
directly while the `preconditioner` will be used to solve the interface problem
via preconditioned conjugate projected gradient algorithm.
"""
function ClassicalFETI(
    mesh, assemble::T; preconditioner::Symbol = :trace) where T <: Base.Callable
    # Validate inputs 
    valid_preconditioners = [:dirichlet, :lumped, :trace]
    if preconditioner ∉  valid_preconditioners
        throw("preconditioner must be in $(valid_preconditioners)")
    end

    local_problems = LocalProblems(mesh, assemble)

    interface_preconditioner = if preconditioner == :trace
        TracePreconditioner() 
    elseif preconditioner == :dirichlet
        DirichletPreconditioner()
    elseif preconditioner == :lumped
        LumpedPreconditioner()
    end

    interface_solver = PCPG()    

    new{typeof(local_problems), 
        typeof(connectivity), 
        typeof(interface_preconditioner),
        typeof(interface_solver)}(; 
        local_problems, connectivity, interface_preconditioner, interface_solver)
end

"""
    solve(...)

Solve a domain decomposed linear static problem using the [`ClassicalFETI`](@ref)
method.

TODO!
"""
function solve(feti::ClassicalFETI, niters, abstol = 1e-6; kwargs...)
    converged = false  
    while !converged && k < niters
        # direct solves on subdomains

        # interface solve 
        
        # Convergence check
        # check_convergend && residual <= abstol && (converged = true)
    end 
end 
