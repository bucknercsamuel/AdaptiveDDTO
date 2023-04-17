#= Adaptive-DDTO -- Parameter structures & functions.
Author: Samuel Buckner (UW-ACL)
=#

# ..:: Quadcopter Object ::..

"""
`Lander` holds the quadcopter parameters.
"""
mutable struct Lander

    # >> Environmental parameters <<
    g::CVector     # [m/s²] Acceleration due to gravity
    ρ::CReal       # [kg/m^3] Air density

    # >> Vehicle parameters <<
    n_rotor::Int   # Number of quadcopter rotors
    mass::CReal    # [kg] Mass of vehicle
    ρ_min::CReal   # [N] Minimum thrust
    ρ_max::CReal   # [N] Maximum thrust

    # >> Constraint parameters <<
    γ_gs::CReal    # [rad] Maximum approach angle
    γ_p::CReal     # [rad] Maximum pointing angle
    v_max_V::CReal # [m/s] Maximum vertical velocity
    v_max_L::CReal # [m/s] Maximum lateral velocity
    τ_max::CReal   # [s] Maximum deferrability between branch points (in terms of node count)

    # >> Initial condition state <<
    r0::CVector    # [m] Initial position
    v0::CVector    # [m] Initial velocity

    # >> Target conditions <<
    n_targs::Int          # Current number of targets
    n_targs_min::Int      # Minimum number of targets
    n_targs_max::Int      # Maximum number of targets
    R_targs::CVector      # [m] Current bounding radii of all targets
    R_targs_min::CReal    # [m] Minimum safe bounding radius for a target
    rf_targs::CMatrix     # [m] Position of each target
    vf_targs::CMatrix     # [m] Velocity of each target
    N_targs::Vector{Int}  # Horizon length of each target
    λ_targs::Vector{Int}  # Order of target rejection
    T_targs::Vector{Int}  # Tag for each target
    ϵ_targs::CVector      # Optimality tolerance
    p_targs::Dict         # Hyperparameters for each target (pcd, prox_veh, prox_clust, µ_99)
    w_des::CVector        # Desirability score weights (pcd, prox_veh, prox_clust, µ_99, R_targs)

    # >> Dynamics <<
    Δt::CReal      # Discretization time-step
    n::Int         # Number of states
    m::Int         # Number of inputs
    A_c::CMatrix   # Continuous-time dynamics A matrix
    B_c::CMatrix   # Continuous-time dynamics B matrix
    p_c::CVector   # Continuous-time dynamics p vector
end

"""
    Lander()
Constructor for the quadcopter.
# Returns
- `lander`: the quadcopter definition.
"""
function Lander()::Lander

    # >> Environmental parameters <<
    g = -9.807*e_z
    ρ = 1.225

    # Rotor Parameters (not stored)
    # (See line 21 of https://github.com/microsoft/AirSim/blob/master/AirLib/include/vehicles/multirotor/RotorParams.hpp)
    C_T = 0.109919
    RPM_max = 6396.667
    d_prop = 0.2286

    # >> Vehicle parameters <<
    n_rotor = 4
    mass = 1
    T_max = n_rotor * C_T * ρ * (RPM_max/60)^2 * d_prop^4 # [N] Max physical thrust of single engine
    ρ_min = 0.2 * T_max # 20% throttle
    ρ_max = 1.0 * T_max # 100% throttle

    # >> Constraint parameters <<
    γ_gs = 89 * DEG_2_RAD
    γ_p = 5 * DEG_2_RAD
    v_max_V = 5
    v_max_L = 5
    τ_max = 5
    # τ_max = 10000 # arbitrarily-large number to disable deferrability constraint

    # >> Initial condition state <<
    r0 = 0*e_x + 0*e_y + 100*e_z
    v0 = 0*e_x + 0*e_y + 0*e_z

    # >> Target conditions <<
    n_targs_min = 3
    n_targs_max = 7
    R_targs_min = 1.
    n_targs = 0 # Need to acquire targets first
    R_targs = CVector(undef, n_targs)
    rf_targs = CMatrix(undef, 3, n_targs)
    vf_targs = zeros(3, n_targs)
    N_targs = Vector{Int}(undef, n_targs)
    λ_targs = Vector{Int}(undef, n_targs)
    T_targs = 1:n_targs
    ϵ_targs = CVector(undef, n_targs)

    # >> Other hyperparameters for each target
    p_targs = Dict(
        "pcd" => CVector(undef, n_targs),         # Point cloud density
        "prox_veh" => CVector(undef, n_targs),    # Proximity of landing site to vehicle
        "prox_clust" => CVector(undef, n_targs),  # Proximity to other landing sites ("cluster proximity")
        "µ_99" => CVector(undef, n_targs),        # 99th percentile uncertainty
    )

    # >> Desirability score weighting <<
    w_des = [0,0,1,0,0] # Weights for: [pcd, prox_veh, prox_clust, µ_99, R_targs]

    # >> Dynamics <<
    Δt = 0.5
    A_c = CMatrix([
        zeros(3,3) I(3);
        zeros(3,3) zeros(3,3)
    ])
    B_c = CMatrix([
        zeros(3,3) zeros(3);
        I(3)/mass  zeros(3)
    ])
    p_c = CVector(vcat(zeros(3),g))
    n,m = size(B_c)

    # >> Make quadcopter object <<
    lander = Lander(
        g,
        ρ,
        n_rotor,
        mass,
        ρ_min,
        ρ_max,
        γ_gs,
        γ_p,
        v_max_V,
        v_max_L,
        τ_max,
        r0,
        v0,
        n_targs,
        n_targs_min,
        n_targs_max,
        R_targs,
        R_targs_min,
        rf_targs,
        vf_targs,
        N_targs,
        λ_targs,
        T_targs,
        ϵ_targs,
        p_targs,
        w_des,
        Δt,
        n,
        m,
        A_c,
        B_c,
        p_c,
    )

    return lander
end


# ..:: Solver Solutions ::..

"""
`Solution` stores the optimal solution.
"""
mutable struct Solution
    # >> Raw data <<
    t::CVector      # [s] Time vector
    r::CMatrix      # [m] Position trajectory
    v::CMatrix      # [m/s] Velocity trajectory
    T::CMatrix      # [m/s^2] Thrust vector
    Γ::CVector      # [m/s^2] Slack thrust magnitude
    cost::CReal     # Optimization's optimal cost

    # >> Processed data <<
    T_nrm::CVector  # [N] Thrust magnitude
    γ::CVector      # [rad] Pointing angle
end

"""
    FailedSolution()
Constructor for a failure solution.
# Arguments
- `sol`: a standard "failed" solution.
"""
function FailedSolution()::Solution

    t = CVector(undef,0)
    r = CMatrix(undef,0,0)
    v = CMatrix(undef,0,0)
    T = CMatrix(undef,0,0)
    Γ = CVector(undef,0)
    cost = Inf
    T_nrm = CVector(undef,0)
    γ = CVector(undef,0)

    return Solution(t,r,v,T,Γ,cost,T_nrm,γ)
end

"""
`DDTOSolution` stores the optimal solution from DDTO.
"""
mutable struct DDTOSolution
    targ_sols::Vector{Solution} # Contains the `Solution` to each target
    costs_sol::CVector          # Costs for each target
    cost_dd::CReal              # Cost for deferred decision
    idx_dd::Int                 # Deferred decision branch point index
end

"""
    EmptyDDTOSolution()
Constructor for an empty DDTOSolution.
# Arguments
- :in n_targs: Number of targets for the solution
- :out sol: a standard empty solution.
"""
function EmptyDDTOSolution(n_targs)::DDTOSolution

    targ_sols = Vector{Solution}(undef, n_targs)
    costs_sol = CVector(undef, n_targs)
    cost_dd   = 0
    idx_dd    = 0

    # Initialize each `Solution` with `FailedSolution()`
    for j = 1:n_targs
        targ_sols[j] = FailedSolution()
    end

    return DDTOSolution(targ_sols,costs_sol,cost_dd,idx_dd)
end

"""
`BranchSolution` stores the solution of a branch from DDTO
"""
mutable struct BranchSolution
    sol::Solution  # Contains the `Solution` for the branch
    cost_dd::CReal # Cost for deferred decision
    idx_dd::Int    # Deferred decision branch point index
    k_comp::Int    # DDTO computation number corresponding to this branch
end