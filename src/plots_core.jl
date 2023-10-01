#= Adaptive-DDTO -- Plotting core functions.
Author: Samuel Buckner (UW-ACL)
=#

function plot_addto_parametric_3D_trajectory(
    lander::Lander, 
    sim::Solution, 
    guid_branches::Vector{Vector{BranchSolution}}, 
    guid_trajs::Array{Solution}; 
    view_az=traj3d_azim)
    """
    Plot a 3D parametric trajectory of the landing vehicle.

    Args:
        lander (Lander): Lander object
        sim (Solution): Simulation solution
        guid_branches (Vector{Vector{BranchSolution}}): Guidance branches
        guid_trajs (Array{Solution}): Guidance trajectories
        view_az (Float64): Azimuthal view angle
    """

    # ..:: Setup ::..
    # Create figure
    fun_name = nameof(var"#self#")
    fig = plt.figure(facecolor="white", figsize=[8,8])
    plt.clf()
    subplot_cols = 1

    # Create and format subplot
    ax = fig.add_subplot(projection="3d", computed_zorder=false)
    ax.set_facecolor("white")
    ax.grid(false)
    ax.xaxis.pane.set_edgecolor("lightgray")
    ax.yaxis.pane.set_edgecolor("lightgray")
    ax.zaxis.pane.set_edgecolor("lightgray")
    ax.xaxis.pane.set_facecolor("white")
    ax.yaxis.pane.set_facecolor("white")
    ax.zaxis.pane.set_facecolor("lightgray")
    ax.view_init(elev=traj3d_elev, azim=view_az)
    ax.dist = 11

    # ..:: Plot guidance trajectory segments ::..
    n_traj_comps = length(guid_trajs)
    ddto_colors = colormap_targs(n_traj_comps, colorant"red", colorant"yellow")
    for k = 1:n_traj_comps
        color = ddto_colors[k]
        ax.plot(guid_trajs[k].r[1,:], guid_trajs[k].r[2,:], guid_trajs[k].r[3,:]; color=color, linewidth=6, alpha=0.5, label="DDTO Sol "*string(k), zorder=10)
    end

    # ..:: Plot DDTO branches ::..
    color = "black"
    branch_iter = 1
    for branch_set in guid_branches
        for branch in branch_set
            color = ddto_colors[branch_iter]
            idx_dd = branch.idx_dd
            r_branch = branch.sol.r
            idx_lb = max(idx_dd, 1)
            idx_ub = idx_lb + min(length(r_branch[1,idx_lb+1:end]), ddto_branch_truncation_idx)
            ax.plot(r_branch[1,idx_lb:idx_ub], r_branch[2,idx_lb:idx_ub], r_branch[3,idx_lb:idx_ub]; color=color, linewidth=1.5, zorder=10)

            # Add truncation indicator
            if idx_ub < length(r_branch[1,:])-1
                ax.plot(r_branch[1,idx_ub], r_branch[2,idx_ub], r_branch[3,idx_ub]; color="red", linestyle="none", marker="s", markersize=3, markeredgecolor="none", zorder=10)
            end
        end
        branch_iter += 1
    end

    # ..:: Plot simulation solution ::..
    ax.plot(sim.r[1,:], sim.r[2,:], sim.r[3,:]; style_ct..., label="Sim", zorder=10)

    # Extra formatting
    ax.set_xlabel("Downrange [m]")
    ax.set_ylabel("Crossrange [m]")
    ax.set_zlabel("Altitude [m]")
    ax.legend(ncol=subplot_cols, loc="center left")
    ax.margins(z=0)
    set_axes_equal(ax)
    plt.tight_layout()

    # ..:: Plot shadows for simulated trajectory and branches ::..
    xlims = ax.get_xlim3d()
    ylims = ax.get_ylim3d()
    x_shadow_offset = xlims[1] # Needs to be adjusted based on viewing azimuth unfortunately
    y_shadow_offset = ylims[2] # Needs to be adjusted based on viewing azimuth unfortunately

    color = "black"
    branch_iter = 1
    for branch_set in guid_branches
        for branch in branch_set
            color = ddto_colors[branch_iter]
            idx_dd = branch.idx_dd
            r_branch = branch.sol.r
            idx_lb = max(idx_dd, 1)
            idx_ub = idx_lb + min(length(r_branch[1,idx_lb+1:end]), ddto_branch_truncation_idx)
            ax.plot(fill(x_shadow_offset, length(r_branch[2,idx_lb:idx_ub])), r_branch[2,idx_lb:idx_ub], r_branch[3,idx_lb:idx_ub]; color="silver", linewidth=1.5, zorder=1)
            ax.plot(r_branch[1,idx_lb:idx_ub], fill(y_shadow_offset, length(r_branch[2,idx_lb:idx_ub])), r_branch[3,idx_lb:idx_ub]; color="silver", linewidth=1.5, zorder=1)
        end
        branch_iter += 1
    end
    ax.plot(fill(x_shadow_offset, length(sim.t)), sim.r[2,:], sim.r[3,:]; style_ct..., color="gray", zorder=1)
    ax.plot(sim.r[1,:], fill(y_shadow_offset, length(sim.t)), sim.r[3,:]; style_ct..., color="gray", zorder=1)

    # Save and show figure
    fig.savefig("$fig_path/$fun_name.$file_ext")
    ;

end

function plot_addto_signals_subplot(
    lander::Lander, 
    sim::Solution, 
    target_radii_history::CMatrix, 
    guid_update_times::CVector,
    guid_lock_time::CReal
)
    """
    Plot various ADDTO signals over time.

    Arguments:
        lander: Lander object
        sim: Simulation solution
        target_radii_history: History of target radii
        guid_update_times: ADDTO guidance update times
        guid_lock_time: ADDTO guidance lock time
    """

    # ..:: Setup ::..
    # Create figure
    fun_name = nameof(var"#self#")
    fig = plt.figure(facecolor="white", figsize=(13,5))
    plt.clf()

    # ..:: Subplot: Thrust Magnitude ::..
    # Create and format subplot
    ax = plt.subplot2grid(shape=(2,2), loc=(0,0), facecolor="white")
    ax.set_facecolor("white")
    ax.grid()

    # Plot simulation solution
    T_nrm_sim = CVector([norm(sim.T[:,i]) for i=1:length(sim.T[1,:])])
    ax.plot(sim.t, T_nrm_sim; style_ct..., label="Sim")

    # Add constraint bounds
    ax.plot(sim.t, fill(lander.ρ_max, length(sim.t)); style_constraint...)
    ax.fill_between(sim.t, lander.ρ_max, lander.ρ_max+3; style_constraint_fill...)
    ax.plot(sim.t, fill(lander.ρ_min, length(sim.t)); style_constraint...)
    ax.fill_between(sim.t, 0, lander.ρ_min; style_constraint_fill...)

    # Extra formatting
    # ax.set_xlabel("Time [s]")
    ax.set_ylabel("Thrust [N]")
    ax.margins(y=0)
    ax.set_xlim([0, sim.t[end]])


    # ..:: Subplot: Tilt Angle ::..
    # Create and format subplot
    ax = plt.subplot2grid(shape=(2,2), loc=(1,0), facecolor="white")
    ax.set_facecolor("white")
    ax.grid()

    # Plot simulation solution
    RAD2DEG = 180/3.14159
    ax.plot(sim.t, sim.γ.*RAD2DEG; style_ct..., label="Sim")

    # Add constraint bounds
    ax.plot(sim.t, fill(lander.γ_p.*RAD2DEG, length(sim.t)); style_constraint...)
    ax.fill_between(sim.t, lander.γ_p.*RAD2DEG, lander.γ_p.*RAD2DEG+1; style_constraint_fill...)

    # Extra formatting
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Tilt Angle ["*L"^\circ"*"]")
    ax.margins(y=0)
    ax.set_xlim([0, sim.t[end]])


    # ..:: Subplot: Radii History ::..
    subplot_cols = 1

    # Create and format subplot
    # ax = fig.add_subplot(1, 2, 2, facecolor="white")
    ax = plt.subplot2grid(shape=(2,2), loc=(0,1), facecolor="white", rowspan=2)
    ax.set_facecolor("white")
    ax.grid()

    # Plotting
    # Radii plots
    colors = colormap_targs(lander.n_targs_max, colorant"magenta", colorant"cyan")
    for j = 1:lander.n_targs_max
        ax.plot(sim.t, target_radii_history[j,:]; style_ct..., color=colors[j], linewidth=3, label="Target "*string(j))
    end

    # DDTO recomputation indicators
    for k = 1:length(guid_update_times)
        if k == 1
            label = "DDTO Update"
        else
            label = ""
        end
        ax.axvline(guid_update_times[k]; color="black", linestyle="dashed", label=label)
    end
    ax.axvline(guid_lock_time; color="limegreen", linestyle="dashed", label="Guid. Lock")

    # Add lower bound constraint
    ax.plot(sim.t, fill(lander.R_targs_min, length(sim.t)); style_constraint...)
    ax.fill_between(sim.t, 0, lander.R_targs_min; style_constraint_fill...)

    # Extra formatting
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Target Radii [m]")
    ax.margins(x=0, y=0)
    ax.legend(ncol=subplot_cols, loc="lower right")
    fig.tight_layout()

    R_max = 0
    for i = 1:size(target_radii_history, 1)
        for j = 1:size(target_radii_history, 2)
            if target_radii_history[i,j] != Inf && target_radii_history[i,j] > R_max
                R_max = target_radii_history[i,j]
            end
        end
    end
    ax.set_ylim([0, R_max*1.1])


    # ..:: Extra Stuff ::..
    # Save and show figure
    fig.savefig("$fig_path/$fun_name.$file_ext")
    ;
end

function replay_addto_landing(
    lander::Lander, 
    sim::Solution, 
    guid_branches::Vector{Vector{BranchSolution}}, 
    guid_trajs::Array{Solution},
    guid_update_times::CVector,
    targs_status::CMatrix,
    targs_radii::CMatrix,; 
    view_az=traj3d_azim)
    """
    Animates an ADDTO landing sequence.

    Arguments:
        lander::Lander: Lander object
        sim::Solution: Simulation solution
        guid_branches::Vector{Vector{BranchSolution}}: Guid. branch solutions
        guid_trajs::Array{Solution}: Guid. trajectories
        guid_update_times::CVector: Guid. update times
        targs_status::CMatrix: Target status history
        targs_radii::CMatrix: Target radii history
        view_az::Float64: Azimuthal view angle for 3D plot
    """

    # ..:: Plot setup ::..
    # Create figure
    fun_name = nameof(var"#self#")
    fig = plt.figure(facecolor="white", figsize=[6,6])
    plt.clf()
    subplot_cols = 1

    # Create and format subplot
    ax = fig.add_subplot(projection="3d", computed_zorder=false)
    ax.set_facecolor("white")
    ax.grid(false)
    ax.xaxis.pane.set_edgecolor("lightgray")
    ax.yaxis.pane.set_edgecolor("lightgray")
    ax.zaxis.pane.set_edgecolor("lightgray")
    ax.xaxis.pane.set_facecolor("white")
    ax.yaxis.pane.set_facecolor("white")
    ax.zaxis.pane.set_facecolor("lightgray")
    ax.view_init(elev=traj3d_elev, azim=view_az)
    ax.dist = 11

    # Extra formatting
    ax.set_title("Adaptive-DDTO: Simulated Landing")
    ax.set_xlabel("Downrange [m]", labelpad=10)
    ax.set_ylabel("Crossrange [m]", labelpad=10)
    ax.set_zlabel("Altitude [m]", labelpad=10)
    ax.margins(z=0)
    ax.set_aspect("equal")

    # Set limits to be equal based off starting altitude (assumes strict descent)
    h0 = sim.r[3,1]
    span = 1.1*h0
    ax.set_xlim([-span/2, span/2])
    ax.set_ylim([-span/2, span/2])
    ax.set_zlim([0, span])

    # Add missing pane frame edge (bug with matplotlib)
    # (uncomment to see what this adds)
    function lims(mplotlims)
        scale = 1.021
        offset = (mplotlims[2] - mplotlims[1])*scale
        return mplotlims[2] - offset, mplotlims[1] + offset
    end
    xlims, ylims, zlims = ax.get_xlim(), ax.get_ylim(), ax.get_zlim()
    xlims_, ylims_, zlims_ = lims(xlims), lims(ylims), lims(zlims)
    iii = transpose([xlims_[1], ylims_[1], zlims_[1]])
    fff = transpose([xlims_[1], ylims_[1], zlims_[2]])
    p = mpl.mplot3d.art3d.Poly3DCollection([[iii, fff]])
    p.set_color("lightgray")
    ax.add_collection3d(p)

    # ..:: Animation design parameters ::..
    alpha_decay = 0.8 # fading out rate for DDTO trajectory spawns
    ds_factor   = 60 # downsampling factor 
    # end_padding = 3 # [s] number of seconds to hold after finishing animation
    time_mult   = 1 # Animation time scaling factor (realtime = 1)
    quad_width  = 15. # [m] width of quadcopter body
    ddto_colors = colormap_targs(lander.n_targs_max, colorant"magenta", colorant"cyan")

    # Global variable initialization (for dynamically-updated variables in animate function)
    global guid_update_iters = 0
    global previous_targ_status = targs_status[:,1]
    global targ_fade_branch = Vector{BranchSolution}(undef, lander.n_targs_max)
    global targ_fade_radii = zeros(lander.n_targs_max)
    global targ_fade_alpha = ones(lander.n_targs_max)
    global next_guid_update_time = guid_update_times[guid_update_iters+1]

    # Variable initialization
    N_sim = length(sim.t)
    dt_sim = sim.t[2] - sim.t[1]
    speedup_ratio = Int(round((1/dt_sim)/ds_factor))


    # ..:: Dynamically-updated plot handles ::..
    global timer_text = ax.text2D(0.5, -0.05, "Sim Time = " * string(0) * " sec", transform=ax.transAxes, fontsize=14, verticalalignment="bottom", horizontalalignment="center")
    global ddto_trajs = [ax.plot([],[],[], color=ddto_colors[k], linewidth=1.5, alpha=1)[1] for k = 1:lander.n_targs_max]
    global ddto_trajs_fading = [ax.plot([],[],[], color="red", linewidth=1.5, alpha=1)[1] for k = 1:lander.n_targs_max] # Used for fade-out effect
    global guid_traj = ax.plot([],[],[], color="gray", linestyle="-", linewidth=2)[1]
    global sim_traj = ax.plot([],[],[], color="black", linewidth=2)[1]
    global ddto_branch_points = [ax.plot([],[],[], color="none", markerfacecolor=ddto_colors[k], markeredgecolor="gray", markeredgewidth=0.5, marker="d", markersize=4)[1] for k = 1:lander.n_targs_max]
    global ddto_branch_points_fading = [ax.plot([],[],[], color="none", markerfacecolor="red", markeredgecolor="none", markeredgewidth=0.5, marker="d", markersize=4)[1] for k = 1:lander.n_targs_max] # Used for fade-out effect
    global quad_body_arm1 = ax.plot([],[],[], color="red", linewidth=1)[1]
    global quad_body_arm2 = ax.plot([],[],[], color="red", linewidth=1)[1]
    global quad_body_center = ax.plot([],[],[], color="none", markerfacecolor="maroon", markeredgecolor="none", marker=".", markersize=7)[1]
    quad_body_legend = ax.plot([],[],[], color="red", markerfacecolor="maroon", markeredgecolor="none", marker=".", markersize=10)[1]

    # Create legend
    midpoint = Int(ceil(lander.n_targs_max/2))
    legend_handles = [
        quad_body_legend,
        sim_traj,
        guid_traj,
        (ddto_trajs[1], ddto_trajs[midpoint], ddto_trajs[lander.n_targs_max]),
        (ddto_branch_points[1], ddto_branch_points[midpoint], ddto_branch_points[lander.n_targs_max]),
    ]
    legend_labels = [
        "Quadcopter",
        "Executed Trajectory",
        "Deferred Trajectory",
        "DDTO Branches",
        "DDTO Branch Points",
    ]

    # Need to use PyCall because handler_map needs Python's native tuple handle :/
    py"""
    import matplotlib
    def custom_legend(ax, handles, labels):
        ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1), ncol=2, frameon=False, handler_map={tuple: matplotlib.legend_handler.HandlerTuple(None)})
    """
    py"custom_legend"(ax, legend_handles, legend_labels)

    # ..:: Define the quadcopter visual geometry at a given point ::..
    function quad_body_geom(x::CReal, y::CReal, z::CReal, width::CReal)::Tuple{CMatrix,CMatrix}
        arm1_x = [x, x]
        arm1_y = [y-width/2, y+width/2]
        arm1_z = z*ones(2)
        arm1 = transpose(hcat(arm1_x, arm1_y, arm1_z))

        arm2_x = [x-width/2, x+width/2]
        arm2_y = [y, y]
        arm2_z = z*ones(2)
        arm2 = transpose(hcat(arm2_x, arm2_y, arm2_z))

        return (arm1,arm2)
    end

    # ..:: Define the animation function ::..
    function animate(iter::Int)
    
        global previous_targ_status, targ_fade_alpha, targ_fade_branch, guid_update_iters, next_guid_update_time
    
        idx = iter * speedup_ratio + 1
    
        if idx <= N_sim-1
            
            # Update timer text
            cur_time = t_sim[idx]
            timer_text.set_text("Sim Time = " * string(round(cur_time, digits=2)) * " sec")

            # Update quad body
            (arm1,arm2) = quad_body_geom(sim.r[1,idx], sim.r[2,idx], sim.r[3,idx], quad_width)
            quad_body_center.set_data([sim.r[1,idx]],[sim.r[2,idx]])
            quad_body_center.set_3d_properties([sim.r[3,idx]])
            quad_body_arm1.set_data(arm1[1,:], arm1[2,:])
            quad_body_arm1.set_3d_properties(arm1[3,:])
            quad_body_arm2.set_data(arm2[1,:], arm2[2,:])
            quad_body_arm2.set_3d_properties(arm2[3,:])

            # Update sim trajectory up to this point
            sim_traj.set_data(sim.r[1,1:idx], sim.r[2,1:idx])
            sim_traj.set_3d_properties(sim.r[3,1:idx])

            # Display guidance lock if at the last guid update
            if cur_time >= guid_update_times[end]
                for k = 1:lander.n_targs_max
                    ddto_trajs[k].set_data([],[])
                    ddto_trajs[k].set_3d_properties([])
                end
                guid_traj.set_data(guid_trajs[end].r[1,:], guid_trajs[end].r[2,:])
                guid_traj.set_3d_properties(guid_trajs[end].r[3,:])

            # Update DDTO trajectories if there is a guid update before the last
            elseif cur_time >= next_guid_update_time      
                # Update new DDTO trajectories
                if guid_update_iters < length(guid_branches)
                    for k = 1:lander.n_targs_max
                        idx_dd = guid_branches[guid_update_iters+1][k].idx_dd
                        ddto_trajs[k].set_data(guid_branches[guid_update_iters+1][k].sol.r[1,idx_dd+1:end], guid_branches[guid_update_iters+1][k].sol.r[2,idx_dd+1:end])
                        ddto_trajs[k].set_3d_properties(guid_branches[guid_update_iters+1][k].sol.r[3,idx_dd+1:end])
                        ddto_branch_points[k].set_data([guid_branches[guid_update_iters+1][k].sol.r[1,idx_dd+1]], [guid_branches[guid_update_iters+1][k].sol.r[2,idx_dd+1]])
                        ddto_branch_points[k].set_3d_properties([guid_branches[guid_update_iters+1][k].sol.r[3,idx_dd+1]])
                    end
                    guid_traj.set_data(guid_trajs[guid_update_iters+1].r[1,:], guid_trajs[guid_update_iters+1].r[2,:])
                    guid_traj.set_3d_properties(guid_trajs[guid_update_iters+1].r[3,:])
                    guid_update_iters += 1
                    next_guid_update_time = guid_update_times[guid_update_iters+1]
                    previous_targ_status = ones(1:lander.n_targs_max)
                end
            end
    
            # Remove targets that have become infeasible
            cur_targ_status = targs_status[:,idx]
            idx_failed_targs = findall(τ->τ<0, cur_targ_status - previous_targ_status)
            for k = 1:length(idx_failed_targs)
                idx_fail = idx_failed_targs[k]
                idx_dd = guid_branches[guid_update_iters][idx_fail].idx_dd
                ddto_trajs[idx_fail].set_data([],[])
                ddto_trajs[idx_fail].set_3d_properties([])
                ddto_trajs_fading[idx_fail].set_data(guid_branches[guid_update_iters][idx_fail].sol.r[1,idx_dd+1:end], guid_branches[guid_update_iters][idx_fail].sol.r[2,idx_dd+1:end])
                ddto_trajs_fading[idx_fail].set_3d_properties(guid_branches[guid_update_iters][idx_fail].sol.r[3,idx_dd+1:end])
                ddto_branch_points[idx_fail].set_data([],[])
                ddto_branch_points[idx_fail].set_3d_properties([])
                ddto_branch_points_fading[idx_fail].set_data([guid_branches[guid_update_iters][idx_fail].sol.r[1,idx_dd+1]], [guid_branches[guid_update_iters][idx_fail].sol.r[2,idx_dd+1]])
                ddto_branch_points_fading[idx_fail].set_3d_properties([guid_branches[guid_update_iters][idx_fail].sol.r[3,idx_dd+1]])
                targ_fade_branch[idx_fail] = guid_branches[guid_update_iters][idx_fail]
                targ_fade_radii[idx_fail] = targs_radii[idx_fail,idx-speedup_ratio]
                targ_fade_alpha[idx_fail] = 1
            end
            previous_targ_status = cur_targ_status

            # Add target zone patches
            if iter > 0
                for patch in ax.patches
                    patch.remove()
                end
            end
            if cur_time >= guid_update_times[end]
                R_lock_idx = findfirst(τ->τ!=-Inf, targs_radii[:,idx])
                R_lock = targs_radii[R_lock_idx,idx]
                patch = plt.Circle((guid_trajs[end].r[1,end], guid_trajs[end].r[2,end]), R_lock, color="gray", alpha=.8)
                ax.add_patch(patch)
                mpl.mplot3d.art3d.pathpatch_2d_to_3d(patch, z=guid_trajs[end].r[3,end], zdir="z")
            else
                for k = 1:lander.n_targs_max
                    patch = plt.Circle((guid_branches[guid_update_iters][k].sol.r[1,end], guid_branches[guid_update_iters][k].sol.r[2,end]), targs_radii[k,idx], color=ddto_colors[k], alpha=0.5)
                    ax.add_patch(patch)
                    mpl.mplot3d.art3d.pathpatch_2d_to_3d(patch, z=guid_branches[guid_update_iters][k].sol.r[3,end], zdir="z")
                end
            end

            # Apply alpha update for fading out
            for k =1:lander.n_targs_max
                ddto_trajs_fading[k].set_alpha(targ_fade_alpha[k])
                ddto_branch_points_fading[k].set_alpha(targ_fade_alpha[k])
                if isassigned(targ_fade_branch,k) && targ_fade_alpha[k] > 0.1
                    patch = plt.Circle((targ_fade_branch[k].sol.r[1,end], targ_fade_branch[k].sol.r[2,end]), targ_fade_radii[k], color="red", alpha=0.5*targ_fade_alpha[k])
                    ax.add_patch(patch)
                    mpl.mplot3d.art3d.pathpatch_2d_to_3d(patch, z=guid_branches[guid_update_iters][k].sol.r[3,end], zdir="z")
                end
                targ_fade_alpha[k] *= alpha_decay
            end
        end

        # Update statement
        if iter > 0 && iter % 10 == 0
            println("Number of frames processed: $iter")
        end
        return timer_text, quad_body_center, quad_body_arm1, quad_body_arm2, sim_traj, ddto_trajs..., guid_traj
    end

    # ..:: Animate! ::..
    anim = animation.FuncAnimation(
        fig, animate, 
        frames=Int(round(N_sim/speedup_ratio)), 
        interval=Int(round(1000/(ds_factor*time_mult))), 
        blit=true, repeat=false)
    
    # Display to notebook
    # video = anim.to_html5_video()
    # html = IJulia.HTML(video)
    # IJulia.display(html)
    # plt.close()

    # Save as GIF
    anim.save("$fig_path/$fun_name.mp4", writer="ffmpeg", fps=60)
    ;
end