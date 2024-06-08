# Combined connectivity, pcpg, and fem to solve a linear system of equations
# (assumes linear static problem) via FETI
using LinearAlgebra

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

struct LocalProblem{TK, TB}
    Ks::TK
    bs::TB

    nullspace_Ks::TK
    pinv_Ks::TK

    function LocalProblem(Ks, bs)
        nullspace_Ks = nullspace(Ks)
        pinv_Ks = pinv(Ks) 
        new(Ks, bs, nullspace_Ks, pinv_Ks)
    end
end     

stiffness(lp::LocalProblem) = lp.Ks
rhs(lp::LocalProblem) = lp.bs
nullspace(lp::LocalProblem) = lp.nullspace_Ks
pinv(lp::LocalProblem) = lp.pinv_Ks

function ClassicalFETI(mesh, assemble; preconditioner::Symbol = :trace)
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

    # Setup connectivity info 
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

    new(; local_problems, connectivity, interface_preconditioner, interface_solver)
end

