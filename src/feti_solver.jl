# Combined connectivity, pcpg, and fem to solve a linear system of equations
# (assumes linear static problem) via FETI
export local_problems, connectivity, interface_solver, interface_preconditioner

"""
    struct ClassicalFETI

Classical finite element tearing and interconnecting (FETI) workspace.

# References
[1] : Farhat, C. A Method of Finite Element Tearing and 
Interconnecting and its Parallel Solution Algorithm. International Journal for 
Numerical Methods in Engineering. vol 32, 1205-1227. 1991. 
"""
@kwdef struct ClassicalFETI{TLP, TC, TIS, TIP}
    local_problems::TLP
    connectivity::TC
    interface_solver::TIS
    interface_preconditioner::TIP
end

local_problems(feti::ClassicalFETI) = feti.local_problems
connectivity(feti::ClassicalFETI) = feti.connectivity
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

    # Assemble local problems 
    partitions = partition(mesh) 
    nparts = length(partitions)
    part_to_local_problem = Vector{LocalProblem}(undef, nparts)
    for (part_i, part) in enumerate(partitions)
        Ks, bs = assemble(part)
        local_problem = LocalProblem(Ks, bs)
        part_to_local_problem[part_i] = local_problem
    end

    # TODO: This is really the hardest part: Setup connectivity info 
    connectivity = Matrix{Connectivity}(undef, nparts, nparts)
    for j in 1:nparts
        for k in j+1:nparts 
            !isneighbor(mesh, j, k) && continue

            partition_j = partitions[j] 
            partition_k = partitions[k]

            # all dofs
            partition_j_dofs = dofs(partition_j)
            partition_k_dofs = dofs(partition_k)

            # specific dofs 
            interface = intersection(partition_j_dofs, partition_k_dofs)
            interior_j = partition_j_dofs[partition_j_dofs .∉  Ref(interface)]
            interior_k = partition_k_dofs[partition_k_dofs .∉  Ref(interface)]
    
            connectivity[j, k] = Connectivity(interface, interior_j)
            connectivity[k, j] = Connectivity(interface, interior_k)
        end
    end

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
