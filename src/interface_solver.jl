# Contains functions and types for the preconditioned conjugate projected gradient
export PCPG, solve!

abstract type InterfaceSolver end 
initialize!(solver::InterfaceSolver) = throw("notimplemented")
niters(s)::Int = s.niters
converged(s)::Bool = s.converged

"""
    solve!(...)

TODO!
"""
function solve!(solver::TIS) where TIS <: InterfaceSolver 

    # solve
    for k = 1:niters(solver)
        _solve!(solver) 
    end

    return nothing 
end 

"""
    struct PCPG <: InterfaceSolver

Preconditioned conjugate projected gradient (PCPG) for solving the interface
problem. 

# References
[1] : Equation (20) from Farhat, C. A Method of Finite Element Tearing and 
Interconnecting and its Parallel Solution Algorithm. International Journal for 
Numerical Methods in Engineering. vol 32, 1205-1227. 1991. 
"""
struct PCPG{TLP, TP <: Preconditioner, Projection, Tλ, TR} <: InterfaceSolver 
    local_problems::TLP 

    M::TP

    P::Projection

    λ₀::Tλ
    λ::Tλ

    rho::Float64

    r₀::TR
    r::TR 
    s::TR 

    niters::Int
    converged::Bool

    function PCPG(; local_problems, M, niters)
        # Interface force matrix, requires pinv of local problems 
        F_I = nothing 
 
        # Initialize lambda
        λ₀ = nothing  
        λ = λ₀
        
        # Initialize the residual 
        r₀ =  nothing 
        r = r₀
        s = r₀

        # Initialize projection operator/matrix   
        P = nothing

        # initialize rho (dot product of residual vectors)
        rho = 1.0

        converged = false 

        return new{
            typeof(local_problems), typeof(M), typeof(P), typeof(λ), typeof(r)}(
            local_problems, M, P, λ₀, λ, rho, r₀, r, s, niters, converged)
    end
end

preconditioner(solver::PCPG) = solver.M
projection(solve::PCPG) = solver.P
lambda(solver::PCPG) = solver.λ
residual(solver::PCPG) = solver.r
rho(solver::PCPG) = solver.rho
search_direction(solver::PCPG) = solver.s 
apply_interface_force(solver::PCPG, s) = throw("notimplemented")

function initialize!(solver::PCPG)
    solver.λ = λ₀
    solver.r = r₀
    solver.s = r₀ # Farhat1991, eq. (20) 
    converged = false 
    return nothing 
end 

"""
    _solve!(solver::PCPG)

A single iteration for solving the interface problem via the [`PCPG`](@ref).
 
# References
[1] : Equation (20) from Farhat, C. A Method of Finite Element Tearing and 
Interconnecting and its Parallel Solution Algorithm. International Journal for 
Numerical Methods in Engineering. vol 32, 1205-1227. 1991. 

[2] : Julia Code 8.1: CG from Darve, E. and Wootters, M. Numerical Linear 
Algebra with Julia. SIAM. 275-277. 2021.
"""
function _solve!(solver::PCPG)
    M = preconditioner(solver)
    r_k_minus_2 = residual(solver)
    r_k_minus_1 = apply(M, r_k_minus_2) # z = M \ r^{(k-2)} 
    rho_old = rho(solver) 
    rho = dot(r_k_minus_2, r_k_minus_1) 
    solver.rho = rho
    β_k = rho/rho_old 

    s_k_minus_1 = search_direction(solver)
    s_k = r_k_minus_1 + β_k*s_k_minus_1
    P = projection(solver)
    s_k = P*s_k
    solver.s_k = s_k 

    q = apply_interface_force(solver, s_k)
    γ_k = rho/q  

    λ_k_minus_1 = lambda(solver)
    λ = λ_k_minus_1 + γ_k*s_k
    solver.λ = λ 
     
    r_k = r_k_minus_1 - γ_k*q     
    solver.r = r_k

    # TODO: some convergence check here 

    return nothing
end 
