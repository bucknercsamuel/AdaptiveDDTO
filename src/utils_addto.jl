#= Adaptive-DDTO -- Utility functions specific to ADDTO.
Author: Samuel Buckner (UW-ACL)
=#

function solve_ddto_stack(lander::Lander)::Tuple{Array{Solution},Array{DDTOSolution}}
    """
    Solve the whole DDTO problem stack.

    Args:
        lander (Lander): The lander object.

    Returns:
        sols_optimal (Array{Solution}): Vectorized container for all single-target solutions.
        sols_ddto (Array{DDTOSolution}): Vectorized container for all DDTO solutions.
    """

    # ..:: Solve for independently-optimal solutions to each target ::..
    (sols_optimal, N_optimal) = solve_optimal_pdg_all_targets(lander)
    lander.N_targs = N_optimal
    costs_optimal = CVector(zeros(lander.n_targs))
    for k = 1:lander.n_targs
       costs_optimal[k] = sols_optimal[k].cost
    end

    # ..:: Solve for DDTO branching solutions to ALL targets ::..
    sols_ddto = solve_ddto_pdg(lander, sols_optimal, costs_optimal)

    return (sols_optimal,sols_ddto)
end

function extract_trunk_segment(lander::Lander, ddto_sol::Array{DDTOSolution})::Solution
    """
    Extract the trunk (deferrable) segment of a full DDTO solution.

    Args:
        lander (Lander): The lander object.
        ddto_sol (Array{DDTOSolution}): Vectorized container for all DDTO solutions.

    Returns:
        sol_trunk (Solution): Container for trunk (deferrable segment) solution.
    """

    # Initialize trunk containers
    r_trunk = CMatrix(undef, 3, 0)
    v_trunk = CMatrix(undef, 3, 0)
    T_trunk = CMatrix(undef, 3, 0)
    Γ_trunk = CVector(undef, 0)
    cost_trunk = 0

    # Build trunk
    for k = 1:length(ddto_sol)
        idx_dd = ddto_sol[k].idx_dd
        r_trunk = hcat(r_trunk, ddto_sol[k].targ_sols[1].r[:,1:idx_dd])
        v_trunk = hcat(v_trunk, ddto_sol[k].targ_sols[1].v[:,1:idx_dd])
        T_trunk = hcat(T_trunk, ddto_sol[k].targ_sols[1].T[:,1:idx_dd])
        append!(Γ_trunk, ddto_sol[k].targ_sols[1].Γ[1:idx_dd])
        cost_trunk += ddto_sol[k].cost_dd
    end

    # Derived variables
    t_trunk = collect(0:length(Γ_trunk)-1) * lander.Δt
    T_nrm_trunk = CVector([norm(T_trunk[:,i],2) for i=1:length(Γ_trunk)])
    γ_trunk = CVector([acos(dot(T_trunk[:,k],e_z)/norm(T_trunk[:,k],2)) for k=1:length(Γ_trunk)])

    sol = Solution(t_trunk, r_trunk, v_trunk, T_trunk, Γ_trunk, cost_trunk, T_nrm_trunk, γ_trunk)
    return sol
end

function extract_guid_lock_traj(lander::Lander, ddto_sol::Array{DDTOSolution}, λ_defer_idx, T_defer::Int, λ_targs_org::Vector{Int})::Solution
    """
    Extract the guidance-locked segment of a full DDTO solution.

    Args:
        lander (Lander): The lander object.
        ddto_sol (Array{DDTOSolution}): Vectorized container for all DDTO solutions.
        λ_defer_idx (Int): λ (preference order) vector index for the target at which the deferral occurs.
        T_defer (Int): Target tag for the target at which the deferral occurs.
        λ_targs_org (Vector{Int}): Original target preference order when DDTO was originally computed.

    Returns:
        sol_guid (Solution): Container for guidance-locked solution.
    """
    
    # Initialize guidance containers
    r_guid = CMatrix(undef, 3, 0)
    v_guid = CMatrix(undef, 3, 0)
    T_guid = CMatrix(undef, 3, 0)
    Γ_guid = CVector(undef, 0)
    cost_guid = 0

    # Build out guidance (trunk up to deferred branch)
    ddto_sol = deepcopy(ddto_sol)
    for k = 1:λ_defer_idx
        # Go from one branch point to the next
        if k < λ_defer_idx
            idx_truncate = ddto_sol[k].idx_dd
            r_guid = hcat(r_guid, ddto_sol[k].targ_sols[1].r[:,1:idx_truncate])
            v_guid = hcat(v_guid, ddto_sol[k].targ_sols[1].v[:,1:idx_truncate])
            T_guid = hcat(T_guid, ddto_sol[k].targ_sols[1].T[:,1:idx_truncate])
            append!(Γ_guid, ddto_sol[k].targ_sols[1].Γ[1:idx_truncate])
            cost_guid += ddto_sol[k].cost_dd
        # Go all the way to landing once at the deferred branch point
        else
            λ_ddto_branch = λ_targs_org[k:end]
            T_ddto_branch = sort(λ_ddto_branch)
            T_defer_idx = findfirst(T_ddto_branch .== T_defer)
            r_guid = hcat(r_guid, ddto_sol[k].targ_sols[T_defer_idx].r[:,1:end])
            v_guid = hcat(v_guid, ddto_sol[k].targ_sols[T_defer_idx].v[:,1:end])
            T_guid = hcat(T_guid, ddto_sol[k].targ_sols[T_defer_idx].T[:,1:end])
            append!(Γ_guid, ddto_sol[k].targ_sols[T_defer_idx].Γ[1:end])
            cost_guid += ddto_sol[k].targ_sols[T_defer_idx].cost
        end
    end

    # Derived variables
    t_guid = collect(0:length(Γ_guid)) * lander.Δt
    T_nrm_guid = CVector([norm(T_guid[:,i],2) for i=1:length(Γ_guid)])
    γ_guid = CVector([acos(dot(T_guid[:,k],e_z)/norm(T_guid[:,k],2)) for k=1:length(Γ_guid)])

    sol = Solution(t_guid, r_guid, v_guid, T_guid, Γ_guid, cost_guid, T_nrm_guid, γ_guid)
    return sol
end

function remove_ddto_target!(lander::Lander, T_targ::Int)
    """
    Remove a target from the lander parameters.
    * NOTE: This function will modify the lander object.

    Args:
        lander (Lander): The lander object.
        T_targ (Int): Target tag for the target to be removed.
    """

    # Determine indices for removal
    pop_idx_T = findfirst(i->i==T_targ, lander.T_targs)
    pop_idx_λ = findfirst(i->i==T_targ, lander.λ_targs)
    slice_T = collect(1:lander.n_targs)
    slice_λ = collect(1:lander.n_targs)
    deleteat!(slice_T, pop_idx_T)
    deleteat!(slice_λ, pop_idx_λ)

    # Parameter updates for removing the target
    lander.n_targs -= 1
    lander.λ_targs = lander.λ_targs[slice_λ]
    lander.T_targs = lander.T_targs[slice_T]
    lander.N_targs = lander.N_targs[slice_T]
    lander.R_targs = lander.R_targs[slice_T]
    lander.ϵ_targs = lander.ϵ_targs[slice_T]
    lander.rf_targs = lander.rf_targs[:,slice_T]
    lander.vf_targs = lander.vf_targs[:,slice_T]
    for (key,~) in lander.p_targs
        lander.p_targs[key] = lander.p_targs[key][slice_T]
    end
end

function switch_decision(lander::Lander, branch_targ::Int)::Bool
    """
    Determine the switch decision at a branch point.

    Args:
        lander (Lander): The lander object.
        branch_targ (Int): Target tag for the target being considered for deferral at the branch point.

    Returns:
        switch_decision (Bool): True if deferral to branch_targ should take place, false otherwise.
    """

    # Find the other targs that aren't the branch targ
    other_targs = copy(lander.T_targs)
    deleteat!(other_targs, findfirst(i->i==branch_targ, other_targs))

    # Get indices for each target
    other_targ_idx = Array{Int}(undef, length(other_targs))
    for j = 1:length(other_targ_idx)
        other_targ_idx[j] = findfirst(i->i==other_targs[j], lander.T_targs)
    end
    branch_targ_idx = findfirst(i->i==branch_targ, lander.T_targs)

    # Check if desirability score of branch targ is greater than that of all other targs
    des_score = zeros(lander.n_targs)
    for j = 1 : lander.n_targs
        des_score[j] = 
            lander.p_targs["pcd"][j] * lander.w_des[1] + 
            lander.p_targs["prox_veh"][j] * lander.w_des[2] + 
            lander.p_targs["prox_clust"][j] * lander.w_des[3] + 
            lander.p_targs["µ_99"][j] * lander.w_des[4] + 
            lander.R_targs[j] * lander.w_des[5]
    end

    if des_score[branch_targ_idx] > maximum(des_score[other_targ_idx])
        switch = true
    else
        switch = false
    end

    return switch
end