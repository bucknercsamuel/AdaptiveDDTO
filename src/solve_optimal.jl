#= Adaptive-DDTO -- Solver functions for single-target trajectory optimization.
Author: Samuel Buckner (UW-ACL)
=#

function solve_optimal_pdg_all_targets(lander::Lander)::Tuple{Vector{Solution},CVector}
    """
    Solve the optimal landing (PDG) problem for a given lander and all targets.

    Args:
        lander (Lander): The lander object.

    Returns:
        solutions (Vector{Solution}): Vectorized container for all single-target solutions.
        N_optimal (CVector): Vectorized container for optimal time horizon for each target.
    """

    # Define container for each `solve_optimal_pdg_single_target` solution
    solutions = Vector{Solution}(undef, lander.n_targs)
    N_optimal = CVector(undef, lander.n_targs)

    # Obtain solutions for each target
    VERB_OPT && println("\n=== Optimal solutions for each target ===")
    for j = 1:lander.n_targs
        
        # Obtain bisection search bounds
        tf_min = norm(lander.rf_targs[:,j] - lander.r0, 2) / sqrt(lander.v_max_L^2 + lander.v_max_V^2) # Fastest straight-line trajectory
        tf_max = tf_min + 20 # There are no fuel/power constraints to enforce an upper bound, so choose 20 seconds arbitrarily
        N_min  = Int(floor((tf_min)/lander.Δt)) + 1
        N_max  = Int(floor((tf_max)/lander.Δt)) + 1

        # Suboptimality tolerance for bisection search
        tol = 1

        # Solve bisection search for optimal trajectory
        VERB_OPT && @printf("\n=== Bisection Search for Target #%i ===\n", j)
        bisection_fun = (N) -> solve_optimal_pdg_single_target(lander, N, j).cost
        N_opt = bisection_search_min_feasible(bisection_fun, N_min, N_max, tol)
        solutions[j] = solve_optimal_pdg_single_target(lander, N_opt, j)

        # Update lander target horizon length
        N_optimal[j] = N_opt
        
        VERB_OPT && @printf("Target: %i, Cost: %.3f\n", lander.T_targs[j], solutions[j].cost)
    end
    # N_optimal .-= 5

    return (solutions, N_optimal)
end

function solve_optimal_pdg_single_target(lander::Lander, N::Int, j_targ::Int)::Solution
    """
    Solve the optimal landing (PDG) problem for a given lander and single target

    Args:
        lander (Lander): The lander object.
        N (Int): Time horizon (used to obtain the optimal time horizon populated into lander.N_targs)
        j_targ (Int): Target index

    Returns:
        sol (Solution): Container for solution variables
    """

    # ..:: Discrete time interval ::..
    Δt = lander.Δt
    tf = Δt * (N-1)
    t  = CVector(range(0, stop=tf, length=N))
    N_ctrl = N-1 # Number of nodes to apply control constraints for (N-1 for ZOH)


    # ..:: Make the optimization problem ::..

    # >> Optimizer setup <<
    if SOLVER == "ECOS"
        mdl = Model(optimizer_with_attributes(ECOS.Optimizer, "verbose" => 0))
    elseif SOLVER == "MOSEK"
        mdl = Model(Mosek.Optimizer)
        JuMP.set_optimizer_attribute(mdl, "LOG",  0) # disable debugging
        JuMP.set_optimizer_attribute(mdl, "MAX_NUM_WARNINGS", 0) # disable warnings
    else
        error("SOLVER is invalid, please select either ECOS or MOSEK")
    end

    # >> Optimization variables <<
    @variable(mdl, r[1:3,1:N])
    @variable(mdl, v[1:3,1:N])
    @variable(mdl, T[1:3,1:N_ctrl])
    @variable(mdl, Γ[1:N_ctrl])

    # >> Convenience functions <<
    X = (k) -> [r[:,k]; v[:,k]] # State at time index k
    U = (k) -> [T[:,k]; Γ[k]]   # Input at time index k

    # >> Cost function <<
    @objective(mdl, Min, Δt*sum(Γ))


    # ..:: Constraints ::..

    # >> Dynamics <<
    A,B,p = c2d_zoh(lander,Δt)
    @constraint(mdl, [k=1:N-1], X(k+1) .==  A*X(k) + B*U(k) + p)

    # >> Thrust bounds <<
    @constraint(mdl, [k=1:N_ctrl], Γ[k] >= lander.ρ_min)
    @constraint(mdl, [k=1:N_ctrl], Γ[k] <= lander.ρ_max)
    @constraint(mdl, [k=1:N_ctrl], vcat(Γ[k], T[:,k]) in MOI.SecondOrderCone(4))

    # >> Attitude pointing constraint <<
    @constraint(mdl, [k=1:N_ctrl], dot(T[:,k],e_z) >= Γ[k]*cos(lander.γ_p))

    # >> Velocity upper bound <<
    @constraint(mdl, [k=1:N], vcat(lander.v_max_L,v[1:2,k]) in MOI.SecondOrderCone(3))
    @constraint(mdl, [k=1:N], v[3,k] >= -lander.v_max_V)
    @constraint(mdl, [k=1:N], v[3,k] <=  lander.v_max_V)

    # >> Boundary conditions <<
    @constraint(mdl, r[:,1] .== lander.r0)
    @constraint(mdl, v[:,1] .== lander.v0)
    @constraint(mdl, r[:,N] .== lander.rf_targs[:,j_targ])
    @constraint(mdl, v[:,N] .== lander.vf_targs[:,j_targ])


    # ..:: Solve the problem and save the solution ::..

    optimize!(mdl)
    if termination_status(mdl) != MOI.OPTIMAL
        return FailedSolution()
    end

    # Raw data
    r = value.(r)
    v = value.(v)
    T = value.(T)
    Γ = value.(Γ)
    cost = objective_value(mdl)

    # Processed data
    T_nrm = CVector([norm(T[:,i],2) for i=1:N_ctrl])
    γ = CVector([acos(dot(T[:,k],e_z)/norm(T[:,k],2)) for k=1:N_ctrl])

    # LCvx slack tightness check
    ϵ_tight = 1e-3
    ϵ_max = 0
    for k=1:N_ctrl
        ϵ_max = max(ϵ_max, abs(T_nrm[k] - Γ[k]))
    end
    if ϵ_max > ϵ_tight
        println("WARNING: LCvx slack tightness violated [OPT] by $(ϵ_max)!")
    end

    # Package the solution
    sol = Solution(t,r,v,T,Γ,cost,T_nrm,γ)

    return sol
end

function bisection_search_min_feasible(fun::Function, τ_min::Int, τ_max::Int, ϵ_tol::Int)::Int
    """
    Find the minimum feasible solution to a function using bisection search
    to a function (not to be confused with DDTO bisection search,
    which finds the *maximum* feasible solution) 

    Args:
        fun (Function): Function to be evaluated (must take opt variable τ as input and return cost)
        τ_min (Int): Bracket search minimum bound
        τ_max (Int): Bracket search maximum bound
        ϵ_tol (Int): Suboptimality convergence tolerance

    Returns:
        τ (Int): Optimal value of τ
    """

    iter = 1
    while (τ_max - τ_min) > ϵ_tol
        # Update τ
        τ = Int(ceil(0.5*(τ_max + τ_min)))

        # Compute feasible DDTO
        cost = fun(τ)

        # Update τ_max or τ_min based on solution convergence
        if ~isinf(cost)
            τ_max = τ
            solve_status = "Feasible"
        else
            τ_min = τ
            solve_status = "Not Feasible"
        end
        VERB_OPT && @printf("Iteration: %i, τ_min: %i, τ_max: %i -- %s\n", iter, τ_min, τ_max, solve_status)

        # Update iteration count
        iter += 1
    end

    # Set optimal τ
    τ_opt = τ_max + 4

    return τ_opt
end
