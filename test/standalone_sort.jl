# Keep the existing ability to load the PyJulia ranking utilities directly
# after defining Lander, without loading the optimizer or ADDTO branch helpers.
using LinearAlgebra
const CReal = Float64
const CVector = Vector{Float64}
const CMatrix = Matrix{Float64}
const e_x = [1.0, 0.0, 0.0]
const e_y = [0.0, 1.0, 0.0]
const e_z = [0.0, 0.0, 1.0]
const DEG_2_RAD = π / 180

include(joinpath(@__DIR__, "..", "src", "params.jl"))
include(joinpath(@__DIR__, "..", "src", "utils_pyjulia.jl"))

lander = Lander()
lander.n_targs_max = 3
reallocate_targ_dims!(lander)
lander.R_targs .= 2.0
lander.p_targs["pcd"] .= NaN
lander.p_targs["prox_veh"] .= NaN
lander.p_targs["µ_99"] .= NaN
lander.p_targs["prox_clust"] .= [0.9, 0.4, 0.2]
sort_des_score!(lander)
println(join(lander.λ_targs, ","))
