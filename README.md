# FETI

<!---[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jfdev001.github.io/FETI.jl/stable/)--->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jfdev001.github.io/FETI.jl/dev/)
[![Build Status](https://github.com/jfdev001/FETI.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jfdev001/FETI.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

A prototype finite element tearing and interconnecting (FETI) solver. Serial algorithm will be implemented first, and possibly the parallel version will be included. A parallel version with a data oriented approach (e.g., using PartitionedArrays.jl) could be more simple than using MPI.

```
## Pseudocode
# Setup local subdomain problems by using a partitioned mesh 
mesh <- make_mesh()
partitions <- partition_mesh(mesh) 
for partition in partitions: 
    ki * ui = bi <- assemble_local_problem(partition)

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
    
    # All dofs for part i not on the interface are interior dofs    
    all_dofs = get_dofs(partitions[part_i])
    interface_dofs = get_all_interface_dofs(connectivity_info, part_i)
    interior_dofs = all_dofs[all_dofs not in interface_dofs]

# Given the delineation between interior and interface dofs, domain 
# decomposition methods are then concerned with handling the fact that
# local problems are not guaranteed to be non-singular (and thus admit 
# a unique solution)... so one performs manipulations to handle this...
...

# The continuity at the interface is in FETI enforced using langrange
# multipliers and solving a saddle point problem
... 
```

# References

[1] : AMFeti. url: https://github.com/AppliedMechanics/AMfeti

[2] : Farhat, C. A Method of Finite Element Tearing and Interconnecting and its Parallel Solution Algorithm. vol 32, 1205-1227. 1991. 
