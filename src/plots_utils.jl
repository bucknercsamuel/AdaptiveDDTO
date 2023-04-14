#= Adaptive-DDTO -- Plotting utility functions & parameter definitions.

Author: Samuel Buckner (UW-ACL)
=#

# ..:: Plotting Utility Parameters ::..

# Setup
PyPlot.svg(true)
fig_path = "figures"
file_ext = "svg"

# Plot styling dictionaries
style_ct = Dict(:color=>"black",:linewidth=>1.5)
style_dt = Dict(:color=>"limegreen",:linestyle=>"none",:marker=>"o",:markersize=>5,:markeredgecolor=>"none")
style_relax = Dict(:color=>"red",:linestyle=>"none",:marker=>"o",:markersize=>3,:markeredgecolor=>"black")
style_thrust = Dict(:edgecolor=>"red",:length=>10,:normalize=>true,:alpha=>0.4)
style_hover = Dict(:color=>"orange",:linewidth=>1.5)
style_apch_cone = Dict(:color=>"blue",:linewidth=>1, :linestyle=>"dashed")
style_ground = Dict(:facecolor=>"gray",:edgecolor=>"black",:alpha=>0.7)
style_constraint = Dict(:color=>"red",:linestyle=>"--",:linewidth=>1)
style_constraint_fill = Dict(:edgecolor=>"none",:facecolor=>"red",:alpha=>0.25)

# Indexing offset based on interpolation method, 1 for ZOH
ctrl_offset = 1

# 3D plot camera viewing angle parameters
traj3d_elev = 20
traj3d_azim = 70

# Set target colors
colorrange_targs = (n,c1,c2) -> collect(range(c1, stop=c2, length=n))
colormap_targs = (n,c1,c2) -> ["#"*hex(color) for color in colorrange_targs(n,c1,c2)]

# Other parameters
ddto_branch_truncation_idx = 15


# ..:: Utility Functions ::..

function modify_styling_dict(styling_dict::Dict, key::String, value::Any)::Dict
    """
    Modify a parameter of an input styling dict.

    Args:
        styling_dict (Dict): The styling dictionary to modify.
        key (String): The key of the parameter to modify.
        value (Any): The new value of the parameter.
    """
    styling_dict_mod = copy(styling_dict)
    styling_dict_mod[Symbol(key)] = value
    return styling_dict_mod
end

function get_thrust_vec(sol::Solution, k::Int, i::Int, j::Int, scale::CReal=0.4)::CVector
    """
    Compute a thrust vector "arrow" for plotting.

    Args:
        sol (Solution): The solution object.
        k (Int): The discrete time step for which to plot.
        i (Int): The "x" axis.
        j (Int): The "y" axis.
        scale (CReal): Multiplicative scaling factor to apply to the thrust vector.
    """
    T = sol.T[:,k]*scale
    r = sol.r[:,k]
    return [r[i],r[j],r[k],T[i],T[j],T[k]]
end

function set_fonts()::Nothing
    """
    Set the fonts for the plots.
    """
    fig_smaller_sz = 10
    fig_small_sz = 13
    fig_med_sz = 15
    fig_big_sz = 17
    plt.rc("font", size=fig_big_sz, family="serif")
    plt.rc("font", serif="Times New Roman")
    plt.rc("axes", titlesize=fig_big_sz)
    plt.rc("axes", labelsize=fig_big_sz)
    plt.rc("xtick", labelsize=fig_med_sz)
    plt.rc("ytick", labelsize=fig_med_sz)
    plt.rc("grid", alpha=0.4)
    plt.rc("legend", fontsize=fig_smaller_sz)
    plt.rc("figure", titlesize=fig_big_sz)
    plt.rc("figure", dpi=300) 
    plt.rc("figure", figsize = [6.5,5])
    plt.rc("animation", html="html5")
    return nothing
end

function set_axes_equal(ax)
    """
    Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Args:
        ax (matplotlib axis): The axis to modify.
    """
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range  = abs(x_limits[2] - x_limits[1])
    x_middle = mean(x_limits)
    y_range  = abs(y_limits[2] - y_limits[1])
    y_middle = mean(y_limits)
    z_range  = abs(z_limits[2] - z_limits[1])
    z_middle = mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*maximum([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

end