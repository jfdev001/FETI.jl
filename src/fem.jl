# Functions and types for FEM for an example problem (poisson on unit square)
# TODO: GMSH --> METIS --> Ferrite.jl or Gridap.jl

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
