# tests for assembling local problem from stiffness matrix 
# https://ferrite-fem.github.io/Ferrite.jl/stable/examples/stokes-flow/#Geometry-and-mesh-generation-with-Gmsh.jl
# https://ferrite-fem.github.io/Ferrite.jl/stable/examples/heat_equation/

@testset "FEM" begin
    # Load a METIS partitioned simple mesh... or just define boundaries
    # yourself in the GMSH???

    # Define a DOF handler for each partition (adds interpolations)

    # Define the assembly function for each partition 
    
    # Create local problems from mesh and assembly operation 
    lps = LocalProblems(mesh, assemble)
        
    # Lesoinne1998 equation (2) asserts this 
    # Δ = ∑(Bˢuˢ) = 0
    @test sum(connectivity(lps) .* dof_to_global(lps)) == 0
end 

