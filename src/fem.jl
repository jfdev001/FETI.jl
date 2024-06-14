# Functions and types for FEM for assembling local problems (stiffness,  load
# vectors) for each subdomain as well as subdomain connectivity
# TODO: GMSH --> METIS --> Ferrite.jl or Gridap.jl
import LinearAlgebra: pinv, nullspace 
export LocalProblems, local_problems, stiffness, rhs, nullspace, pinv, 
    connectivity, dof_to_global

""" 
    struct LocalProblem 

TODO!
"""
struct LocalProblem{TK, TRHS, TB, TDOF}
    "Connectivity matrix on subdomain `s`"
    connectivity_Bs::TK
    "Stiffness matrix on subdomain `s`"
    stiffness_Ks::TK
    "Right hand side (load vector) on subdomain `s`"
    rhs_bs::TB
    "Nullspace of `Ks` -- also denoted `Rs` in Farhat1991"
    nullspace_Ks::TK
    "Cached pseudoinverse of stiffness matrix `Ks`"
    pinv_Ks::TK
    "Maps a dof in 1:length(rhs_bs) to its global id"
    dof_to_global::TDOF

    function LocalProblem(Ks, bs, Bs)
        nullspace_Ks = nullspace(stiffness_Ks)
        pinv_Ks = pinv(stiffness_Ks) 
        dof_to_global = collect(1:length(bs))  # TODO:
        return new{typeof(Ks), typeof(bs), typeof(Bs), typeof(dof_to_global)}(
            Bs, stiffness_Ks, rhs_bs, nullspace_Ks, pinv_Ks, dof_to_global)
    end
end     

stiffness(lp::LocalProblem) = lp.stiffness_Ks
rhs(lp::LocalProblem) = lp.rhs_bs
nullspace(lp::LocalProblem) = lp.nullspace_Ks
pinv(lp::LocalProblem) = lp.pinv_Ks
connectivity(lp::LocalProblem) = lp.connectivity_Bs
dof_to_global(lp::LocalProblem) = lp.dof_to_global

"""
    struct LocalProblems

TODO!
"""
struct LocalProblems{TLPS}
    local_problems::TLPS

    function LocalProblems(mesh, assemble)
        local_problems = nothing # todo 

        # Assemble local problems 
        global_dofs = dofs(mesh) # should be repeats in here for interfaces 
        partitions = partition(mesh) 
        nparts = length(partitions)
        part_to_local_problem = Vector{LocalProblem}(undef, nparts)
        n_local_dofs = 0

        """
        # Equation (2) from Lesoinne1998 for feti signed matrix    
        # assemble local problems 
        local_probleims = []
        n_I = n_global_interface_dofs = sum(map(dof->hasduplicate(dof), global_dofs))
        for j in 1:nparts # gets all nparts
            part_j = partitions[j]
            Kj, bj = assemble(part_j)
            neighbor_k_to_Bj = Dict{Int8, Matrix{Int8}}()
            # get connectivity... need to figure out shape of this matrix...
            for k in j+1:nparts
                part_k = partitions[k]
                if isneighbor(mesh, j, k) # isneighbor(part_j, part_k)
                    n_s_dofs = ndofs(part)
                    
                end 
            end
            lp = LocalProblem(Kj, bj, Bj)
            push!(local_problems, lp)
        end
        """
        for (part_i, part) in enumerate(partitions) # still need to check all other parts...
            Ks, bs = assemble(part)
            local_to_dof = collect(1:length(bs)) 
            local_problem = LocalProblem(Ks, bs)
            part_to_local_problem[part_i] = local_problem
        end

        # TODO: This is really the hardest part: Setup connectivity info 
        # TODO: The logic below isn't ideal... something more like
        """
        global_dofs::Vector = global_dofs(mesh)
        
        """
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

        return new{typeof(local_problems)}(local_problems) 
    end 
end

local_problems(lps::LocalProblems) = lps.local_problems

# Getters returning fields for all subdomains
stiffness(lps::LocalProblems) = stiffness.(local_problems(lps))
rhs(lps::LocalProblems) = rhs.(local_problems(lps))
nullspace(lps::LocalProblems) = nullspace.(local_problems(lps))
pinv(lps::LocalProblems) = pinv.(local_problems(lps)) 
connectivity(lps::LocalProblems) = connectivity.(local_problems(lps)) 
dof_to_global(lps::LocalProblem) = dof_to_global.(local_problems(lps))

# Getters returning field of a local problem for a subdomain `s`
stiffness(lps::LocalProblems, s::Int) = stiffness.(local_problems(lps))[s]
rhs(lps::LocalProblem, s::Int) = rhs.(local_problems(lps))[s]
nullspace(lps::LocalProblems, s::Int) = nullspace.(local_problems(lps))[s]
pinv(lps::LocalProblems, s::Int) = pinv.(local_problems(lps))[s]
connectivity(lps::LocalProblems, s::Int) = connectivity.(local_problems(lps))[s]
dof_to_global(lps::LocalProblem, s::Int) = dof_to_global.(local_problems(lps))[s]

function partition(mesh)
    throw("notimplemented")
end

function isneighbor(mesh, j, k)::Bool
    throw("notimplemented")
end 

function ndofs(partition)
    throw("notimplemented")
end 

