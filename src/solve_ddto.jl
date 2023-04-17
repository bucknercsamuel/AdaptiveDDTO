#= Adaptive-DDTO -- Solver functions for deferred-decision trajectory optimization (DDTO).
Author: Samuel Buckner (UW-ACL)
=#

function solve_ddto_pdg(lander::Lander, sols_optimal_0::Vector{Solution}, costs_optimal_0::CVector)::Vector{DDTOSolution}
    """
    Top-level DDTO solver for all branch points

    Args:
        lander (Lander): The lander object.
        sols_optimal_0 (Vector{Solution}): Optimal solutions from initial condition.
        costs_optimal_0 (CVector): Optimal costs from initial condition.

    Returns:
        Vector{DDTOSolution}: Vectorized container for all DDTO branch solutions.
    """

    # Define container for each DDTO branch solution
    ddto_branch_sols = Vector{DDTOSolution}(undef, lander.n_targs)
    for k = 1:(lander.n_targs)
        ddto_branch_sols[k] = EmptyDDTOSolution(lander.n_targs-k+1)
    end

    # Define first "previous" ddto solution for first branch using optimal solutions
    previous_ddto_solution = EmptyDDTOSolution(lander.n_targs)
    for k=1:(lander.n_targs)
        previous_ddto_solution.targ_sols[k] = copy(sols_optimal_0[k])
        previous_ddto_solution.costs_sol[k] = copy(costs_optimal_0[k])
    end

    # Define running deferred-decision (DD) trajectory segment cost sum
    cost_dd_sum = 0.

    # Set scheduling to update ϵ_targs (all targets) for each DDTO subproblem (n-1 subproblems in total)
    ϵ_targ_schedule = LinRange(0.1, 1.0, lander.n_targs-1)

    # Solve each DDTO branching subproblem in order of rejection
    n_targs_total = copy(lander.n_targs)
    lander_ = copy(lander) # Temp object to be mutated through DDTO loop
    pop_idx = 0
    for k = 1:(n_targs_total-1)

        # Set the suboptimality tolerances for this iteration
        for j in 1:lander_.n_targs
            lander_.ϵ_targs[j] = ϵ_targ_schedule[k]
        end 

        if VERB_DDTO
            specifiers = repeat("%.3f, ", lander_.n_targs)
            specifiers = specifiers[1:end-2] # Remove string and comma at end
            format_string = "   Chosen suboptimality tolerances: {"*specifiers*"}\n"

            @printf("\n========= Solving DDTO for Branch #%i =========\n", k)
            @eval @printf($format_string, $lander_.ϵ_targs...)
        end

        # Obtain Bisection-optimal DDTO solution for this branch
        if k > 1
            previous_ddto_solution = copy(ddto_branch_sols[k-1])
            previous_ddto_solution.idx_dd = 0
            previous_ddto_solution.cost_dd = 0
        
            # Update previous_ddto_solution to not have the previously-removed target
            deleteat!(previous_ddto_solution.targ_sols, pop_idx)
            deleteat!(previous_ddto_solution.costs_sol, pop_idx)

            # Truncate previous_ddto_solution for previous solution's deferral
            trunc_start = ddto_branch_sols[k-1].idx_dd+1
            for j=1:length(previous_ddto_solution.targ_sols)
                n_nodes = length(previous_ddto_solution.targ_sols[j].t)
                previous_ddto_solution.targ_sols[j].t     = previous_ddto_solution.targ_sols[j].t[1:(n_nodes-trunc_start+1)]
                previous_ddto_solution.targ_sols[j].r     = previous_ddto_solution.targ_sols[j].r[:,trunc_start:end]
                previous_ddto_solution.targ_sols[j].v     = previous_ddto_solution.targ_sols[j].v[:,trunc_start:end]
                previous_ddto_solution.targ_sols[j].T     = previous_ddto_solution.targ_sols[j].T[:,trunc_start:end]
                previous_ddto_solution.targ_sols[j].Γ     = previous_ddto_solution.targ_sols[j].Γ[trunc_start:end]
                previous_ddto_solution.targ_sols[j].γ     = previous_ddto_solution.targ_sols[j].γ[trunc_start:end]
                previous_ddto_solution.targ_sols[j].T_nrm = previous_ddto_solution.targ_sols[j].T_nrm[trunc_start:end]
            end
        end

        # Compute DDTO bisection search
        ddto_branch_sols[k] = solve_bisection_qcvx_ddto(lander_, costs_optimal_0, cost_dd_sum, previous_ddto_solution)
        
        # Determine target to be removed (first in the current list of λ_targs)
        λ_targ = lander_.λ_targs[1]
        deleteat!(lander_.λ_targs, 1)
        pop_idx = findfirst(i->i==λ_targ, lander_.T_targs)

        # Have to do some slicing magic for matrices
        matrix_slice = collect(1:lander_.n_targs)
        deleteat!(matrix_slice, pop_idx)

        # Update lander target and IC properties for removal of deferred target
        lander_.n_targs -= 1
        deleteat!(lander_.T_targs, pop_idx)
        deleteat!(lander_.N_targs, pop_idx)
        deleteat!(lander_.R_targs, pop_idx)
        deleteat!(lander_.ϵ_targs, pop_idx)
        lander_.N_targs .-= ddto_branch_sols[k].idx_dd
        lander_.r0 = ddto_branch_sols[k].targ_sols[1].r[:,ddto_branch_sols[k].idx_dd+1]
        lander_.v0 = ddto_branch_sols[k].targ_sols[1].v[:,ddto_branch_sols[k].idx_dd+1]
        lander_.rf_targs = lander_.rf_targs[:,matrix_slice]
        lander_.vf_targs = lander_.vf_targs[:,matrix_slice]


        # Update deferred-decision (DD) cost for next branch iteration
        cost_dd_sum += ddto_branch_sols[k].cost_dd

        # Parameter update print statements
        if VERB_DDTO && (k < n_targs_total-1)
            @printf("   Removed target %i for next branch iteration\n", λ_targ)
        end
    end

    # Add a final element to the branch solutions for the final target
    ddto_branch_sols[end].targ_sols = [copy(ddto_branch_sols[end-1].targ_sols[2])]
    ddto_branch_sols[end].costs_sol = [copy(ddto_branch_sols[end-1].costs_sol[2])]
    ddto_branch_sols[end].idx_dd    = 0
    ddto_branch_sols[end].cost_dd   = 0

    # The single final solution must be truncated for consistency
    trunc_start = ddto_branch_sols[end-1].idx_dd+1
    n_nodes = length(ddto_branch_sols[end].targ_sols[1].t)
    ddto_branch_sols[end].targ_sols[1].t     = ddto_branch_sols[end].targ_sols[1].t[1:(n_nodes-trunc_start+1)]
    ddto_branch_sols[end].targ_sols[1].r     = ddto_branch_sols[end].targ_sols[1].r[:,trunc_start:end]
    ddto_branch_sols[end].targ_sols[1].v     = ddto_branch_sols[end].targ_sols[1].v[:,trunc_start:end]
    ddto_branch_sols[end].targ_sols[1].T     = ddto_branch_sols[end].targ_sols[1].T[:,trunc_start:end]
    ddto_branch_sols[end].targ_sols[1].Γ     = ddto_branch_sols[end].targ_sols[1].Γ[trunc_start:end]
    ddto_branch_sols[end].targ_sols[1].γ     = ddto_branch_sols[end].targ_sols[1].γ[trunc_start:end]
    ddto_branch_sols[end].targ_sols[1].T_nrm = ddto_branch_sols[end].targ_sols[1].T_nrm[trunc_start:end]

    return ddto_branch_sols
end

function solve_bisection_qcvx_ddto(lander::Lander, costs_optimal::CVector, cost_dd::CReal, previous_ddto_solution::DDTOSolution)::DDTOSolution
    """
    Uses bisection search to solve quasiconvex maximal-deferral optimization problem in DDTO.

    Args:
        lander: The lander object
        costs_optimal: Optimal costs from `solve_optimal_pdg_all_targets`
        cost_dd: Running cost for decision deferral
        previous_ddto_solution: Contains the previous branching subproblem solution

    Returns:
        ddto_solution: Contains the DDTO solution for this target/branch point
    """

    # Initial search bracket
    τ_min = 0
    τ_max = min(min(lander.N_targs...) - 2, lander.τ_max)

    # Bisection search to solve quasiconvex (QCvx) optimization problem
    VERB_DDTO && println("=== Bisection Search for QCvx Optimization ===")
    iter = 1
    while (τ_max - τ_min) > 1
        # Update τ
        τ = Int(ceil(0.5*(τ_max + τ_min)))

        # Compute feasible DDTO
        (~, status_feas) = solve_feasible_ddto(lander, τ, costs_optimal, cost_dd)

        # Update τ_min or τ_max based on solution convergence
        if status_feas == MOI.OPTIMAL
            τ_min = τ
            solve_status = "Feasible"
        else
            τ_max = τ
            solve_status = "Not Feasible"
        end
        VERB_DDTO && @printf("Iteration: %i, τ_min: %i, τ_max: %i -- %s\n", iter, τ_min, τ_max, solve_status)

        # Update iteration count
        iter += 1
    end

    # Set optimal τ
    τ_opt = τ_min
    VERB_DDTO && println("Bisection search terminated -- reached convergence condition (τ_max - τ_min) = 1")

    # Compute converged DDTO solution
    # (just re-use previous solution if cannot be deferred)
    if τ_opt == 0
        ddto_solution = copy(previous_ddto_solution)
        status_feas = MOI.OPTIMAL
    else
        (ddto_solution, status_feas) = solve_feasible_ddto(lander, τ_opt, costs_optimal, cost_dd)
        ddto_solution.idx_dd = τ_opt
    end

    # Determine solution convergence
    if status_feas == MOI.OPTIMAL
        VERB_DDTO && @printf("Bisection search successful -- τ_opt: %i\n", τ_opt)
    elseif status_feas == MOI.SLOW_PROGRESS
        @printf("SOLVER WARNING: Slow progress on solution -- proceeding as if optimal -- τ_opt: %i\n", τ_opt)
    else
        error("SOLVER ERROR: Bisection search unsuccessful. Problem is unsolved.")
    end
    VERB_DDTO && println("Deferred-decision Solutions:")
    for j = 1:lander.n_targs
        VERB_DDTO && @printf("   Target: %i, Cost: %.3f\n", lander.T_targs[j], ddto_solution.costs_sol[j])
    end

    return ddto_solution

end

function solve_feasible_ddto(lander::Lander, τ::Int, costs_optimal::CVector, cost_dd::CReal)::Tuple{DDTOSolution, MOI.TerminationStatusCode}
    """
    Solves the DDTO feasibility problem for a given deferrability index τ in the bisection search.

    Args:
        lander: The lander object
        τ: deferrability index
        costs_optimal: Optimal costs from `solve_optimal_pdg_all_targets`
        cost_dd: Running cost for decision deferral

    Returns:
        ddto_solution: Contains the DDTO solution
        feas_status: Feasibility problem solution status code (see MOI.TerminationStatusCode documentation)

    """

    # ..:: Discrete time interval ::..
    N  = max(lander.N_targs...)
    n  = lander.n_targs
    Δt = lander.Δt
    tf = Δt * (N-1)
    N_ctrl = N-1 # Number of nodes to apply control constraints for (N-1 for ZOH)
    A,B,p = c2d_zoh(lander,Δt)


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
    @variable(mdl, r[1:3,1:N,1:n])
    @variable(mdl, v[1:3,1:N,1:n])
    @variable(mdl, T[1:3,1:N,1:n])
    @variable(mdl, Γ[1:N,1:n])

    # >> Expression holders <<
    subopt = Array{QuadExpr}(undef, N_ctrl, n)

    # >> Convenience functions <<
    X = (k,j) -> [r[:,k,j]; v[:,k,j]] # State at time index k and target j
    U = (k,j) -> [T[:,k,j]; Γ[k,j]]    # Input at time index k and target j

    # >> Cost function <<
    # (Feasibility problem, cost is trivially set to 0)
    @objective(mdl, Min, 0)


    # ..:: Constraints ::..

    # >> Iterate through targets <<
    for j = 1:n

        # Target N
        N_targ = lander.N_targs[j]
        N_targ_ctrl = N_targ - 1

        # Slice indexing to n without current target j
        J = collect(1:n)
        deleteat!(J, j)

        # >> Dynamics <<
        @constraint(mdl, [k=1:N_targ-1], X(k+1,j) .== A*X(k,j) + B*U(k,j) + p)

        # >> Thrust bounds <<
        @constraint(mdl, [k=1:N_targ_ctrl], Γ[k,j] >= lander.ρ_min)
        @constraint(mdl, [k=1:N_targ_ctrl], Γ[k,j] <= lander.ρ_max)
        @constraint(mdl, [k=1:N_targ_ctrl], vcat(Γ[k,j], T[:,k,j]) in MOI.SecondOrderCone(4))

        # >> Attitude pointing constraint <<
        @constraint(mdl, [k=1:N_targ_ctrl], dot(T[:,k,j],e_z) >= Γ[k,j]*cos(lander.γ_p))

        # >> Velocity upper bound <<
        @constraint(mdl, [k=1:N_targ], vcat(lander.v_max_L,v[1:2,k,j]) in MOI.SecondOrderCone(3))
        @constraint(mdl, [k=1:N_targ], v[3,k,j] >= -lander.v_max_V)
        @constraint(mdl, [k=1:N_targ], v[3,k,j] <=  lander.v_max_V)

        # >> Identicality << 
        for k = 1:N_targ_ctrl
            if τ > 0
                if k <= τ
                    for l = 1:n-1
                        @constraint(mdl, U(k,j) .== U(k,J[l]))
                    end
                end
            end
        end

        # >> Suboptimality <<
        for k = 1:N_ctrl
            if k <= N_targ_ctrl
                subopt[k,j] = @expression(mdl, Δt*Γ[k,j])
            else
                subopt[k,j] = @expression(mdl, 0.0)
            end
        end

        # >> Zero out state/control nodes from N_targ+1 to N <<
        @constraint(mdl, [k=N_targ+1:N],           X(k,j) .== zeros(lander.n,1))
        @constraint(mdl, [k=N_targ_ctrl+1:N_ctrl], U(k,j) .== zeros(lander.m,1))

        # >> Boundary conditions << 
        @constraint(mdl, r[:,1,j]      .== lander.r0)
        @constraint(mdl, v[:,1,j]      .== lander.v0)
        @constraint(mdl, r[:,N_targ,j] .== lander.rf_targs[:,j])
        @constraint(mdl, v[:,N_targ,j] .== lander.vf_targs[:,j])

        # >> Sub-optimality <<
        @constraint(mdl, sum(subopt[:,j]) + cost_dd .<= (1 + lander.ϵ_targs[j]) * costs_optimal[j])
    end

    optimize!(mdl)
    feas_status = JuMP.termination_status(mdl)

    r = value.(r)
    v = value.(v)
    T = value.(T)
    Γ = value.(Γ)


    # ..:: Determine optimal cost and deferred-decision (DD) cost ::..

    costs_sol = CVector(zeros(n))
    cost_dd  = 0
    for j = 1:n
        N_targ = lander.N_targs[j]
        for k = 1:N_targ-1
            costs_sol[j] += Δt*Γ[k,j]
            if k==τ && j==1
                cost_dd = costs_sol[j]
            end
        end
    end


    # ..:: Package the DDTO Solution ::..

    ddto_solution = EmptyDDTOSolution(n)
    for j = 1:n
        N_targ = lander.N_targs[j]
        N_targ_ctrl = N_targ - 1

        # Raw data
        ddto_solution.targ_sols[j].t = CVector(range(0, stop=tf, length=lander.N_targs[j]))
        ddto_solution.targ_sols[j].r = r[:,1:N_targ,j]
        ddto_solution.targ_sols[j].v = v[:,1:N_targ,j]
        ddto_solution.targ_sols[j].T = T[:,1:N_targ_ctrl,j]
        ddto_solution.targ_sols[j].Γ = Γ[1:N_targ_ctrl,j]
        ddto_solution.targ_sols[j].cost = costs_sol[j]

        # Processed data
        ddto_solution.targ_sols[j].T_nrm = CVector([norm(T[:,k,j],2) for k=1:N_targ_ctrl])
        ddto_solution.targ_sols[j].γ     = CVector([acos(dot(T[:,k,j],e_z)/norm(T[:,k,j],2)) for k=1:N_targ_ctrl])

        # LCvx slack tightness check
        ϵ_tight = 1e-1
        ϵ_max = 0
        for k=1:N_targ_ctrl
            ϵ_max = max(ϵ_max, abs(ddto_solution.targ_sols[j].T_nrm[k] - Γ[k,j]))
        end
        if ϵ_max > ϵ_tight
            println("WARNING: LCvx slack tightness violated [DDTO] by $(ϵ_max)!")
        end
    end
    ddto_solution.costs_sol = costs_sol
    ddto_solution.cost_dd   = cost_dd

    return (ddto_solution, feas_status)

end