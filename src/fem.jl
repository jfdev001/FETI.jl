# Functions and types for FEM for an example problem (poisson on unit square)
# TODO: GMSH --> METIS --> Ferrite.jl or Gridap.jl
import LinearAlgebra: pinv, nullspace 
export LocalProblem, stiffness, rhs, nullspace, pinv

""" 
    struct LocalProblem 

TODO!
"""
struct LocalProblem{TK, TB}
    "Stiffness matrix on subdomain `s`"
    Ks::TK
    "Right hand side (forcing vector) on subdomain `s`"
    bs::TB
    "Nullspace of `Ks` -- also denoted `Rs` in Farhat1991"
    nullspace_Ks::TK
    "Cached pseudoinverse of stiffness matrix `Ks`"
    pinv_Ks::TK

    function LocalProblem(Ks, bs)
        nullspace_Ks = nullspace(Ks)
        pinv_Ks = pinv(Ks) 
        return new{typeof(Ks), typeof(bs)}(Ks, bs, nullspace_Ks, pinv_Ks)
    end

    function LocalProblem(Ks, bs, nullspace_Ks, pinv_Ks)
        return new{typeof(Ks), typeof(bs)}(Ks, bs, nullspace_Ks, pinv_Ks)
    end
end     

stiffness(lp::LocalProblem) = lp.Ks
rhs(lp::LocalProblem) = lp.bs
nullspace(lp::LocalProblem) = lp.nullspace_Ks
pinv(lp::LocalProblem) = lp.pinv_Ks

"""
    struct Connectivity 

TODO!
"""
struct Connectivity

end 

function partition(mesh)
    throw("notimplemented")
end

function isneighbor(mesh, j, k)::Bool
    throw("notimplemented")
end 

function dofs(partition)
    throw("notimplemented")
end 

