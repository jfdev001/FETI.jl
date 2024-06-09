@testset "FETI Solver" begin
    singular_matrix = [1.5 1.3; 1.2 1.9]      
    b = rand(size(singular_matrix, 2))
    lp1 = FETI.LocalProblem(singular_matrix, b)
    lp2 = FETI.LocalProblem(stiffness(lp1), rhs(lp1), nullspace(lp1), pinv(lp1))
end 
