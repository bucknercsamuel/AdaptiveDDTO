#= Adaptive-DDTO -- Setup script.
Author: Samuel Buckner (UW-ACL)
=#

# Necessary Imports
using LinearAlgebra
using Random
using JuMP, MosekTools, ECOS
using Statistics
using Printf
using PyPlot
using PyCall
using Colors
using Noise
using IJulia

# Data types
const CReal = Float64
const CVector = Vector{CReal}
const CMatrix = Matrix{CReal}

# Inertial frame vectors
const e_x = CVector([1,0,0])
const e_y = CVector([0,1,0])
const e_z = CVector([0,0,1])

# Unit conversion
const RAD_2_DEG = 180 / π
const DEG_2_RAD = π / 180
const M_2_KM = 1/1000
const KM_2_M = 1000
const N_2_KN = 1/1000
const KN_2_N = 1000

# Set solver
SOLVER = "ECOS"

# Set RNG type and randomization seed
rng = MersenneTwister(1234)
Random.seed!(2007)

# Set verbose option for each algorithm
VERB_OPT = false # Choose whether to print internal updates for the optimal solution bracket searches
VERB_DDTO = false # Choose whether to print internal updates for the DDTO solution branches

# Configure copy method to work with Lander structure
Base.copy(t::T) where T = T([deepcopy(getfield(t,k)) for k ∈ fieldnames(T)]...)