#= Adaptive-DDTO -- Utility functions for PyJulia interface (only used for AirSim code).
Author: Samuel Buckner (UW-ACL)
=#

function rk4_step_pyjulia(lander::Lander, guid_traj::Solution, x_cur::CVector, t_cur::CReal, Δt::CReal)::CVector
    """
    Necessary wrapper for rk4_step where the dynamics function is created here
    (since there is no way to pass an anonymous function from Python to Julia AFAIK)

    Args:
        lander (Lander): The lander object.
        guid_traj (Solution): The guidance trajectory.
        x_cur (CVector): Current state.
        t_cur (CReal): Current time.
        Δt (CReal): Time step.

    Returns:
        x_cur (CVector): Updated state.
    """

    # Get current dynamics with controller tracking `guid_traj`
    ct_dyn = (t,x) -> lander.A_c*x + lander.B_c*optimal_controller(t,x,guid_traj) + lander.p_c # Continuous-time dynamics

    # Call `rk4_step`
    x_cur = rk4_step(x_cur, ct_dyn, t_cur, Δt)

    return x_cur
end

function reallocate_targ_dims!(lander::Lander)
    """
    Reallocate the target dimensions to the maximum number of targets.
    * NOTE: This function will modify the lander object.

    Args:
        lander (Lander): the lander object.
    """

    lander.n_targs = lander.n_targs_max
    lander.R_targs = CVector(undef, lander.n_targs)
    lander.rf_targs = CMatrix(undef, 3, lander.n_targs)
    lander.vf_targs = zeros(3, lander.n_targs)
    lander.N_targs = Vector{Int}(undef, lander.n_targs)
    lander.λ_targs = Vector{Int}(undef, lander.n_targs)
    lander.T_targs = 1:lander.n_targs
    lander.ϵ_targs = CVector(undef, lander.n_targs)
    for (key,~) in lander.p_targs
        lander.p_targs[key] = CVector(undef, lander.n_targs)
    end
end

function sort_des_score!(lander::Lander)
    """
    Sort the targets by descending desirability score.
    * NOTE: This function will modify the lander object.

    Args:
        lander (Lander): the lander object.
    """

    des_score = zeros(lander.n_targs)
    for j = 1 : lander.n_targs
        des_score[j] = 
            lander.p_targs["pcd"][j] * lander.w_des[1] + 
            lander.p_targs["prox_veh"][j] * lander.w_des[2] + 
            lander.p_targs["prox_clust"][j] * lander.w_des[3] + 
            lander.p_targs["µ_99"][j] * lander.w_des[4] + 
            lander.R_targs[j] * lander.w_des[5]
    end
    lander.λ_targs = sortperm(des_score)

    # Sort the last two targets to always be increasing since there is no rejection preference between the two
    lander.λ_targs[end-1:end] = sort(lander.λ_targs[end-1:end])
end