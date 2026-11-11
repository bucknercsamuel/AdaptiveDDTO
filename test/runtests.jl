using Test
using LinearAlgebra

# The decision functions need only Julia's standard library; avoid loading the
# optional solver and plotting packages for these regression tests.
const CReal = Float64
const CVector = Vector{Float64}
const CMatrix = Matrix{Float64}
const e_x = [1.0, 0.0, 0.0]
const e_y = [0.0, 1.0, 0.0]
const e_z = [0.0, 0.0, 1.0]
const DEG_2_RAD = π / 180

include(joinpath(@__DIR__, "..", "src", "params.jl"))
include(joinpath(@__DIR__, "..", "src", "utils_addto.jl"))
include(joinpath(@__DIR__, "..", "src", "utils_pyjulia.jl"))

function three_targets(; inactive=NaN)
    lander = Lander()
    lander.n_targs_max = 3
    reallocate_targ_dims!(lander)
    lander.R_targs .= 2.0
    lander.p_targs["pcd"] .= inactive
    lander.p_targs["prox_veh"] .= inactive
    lander.p_targs["µ_99"] .= inactive
    lander.p_targs["prox_clust"] .= [0.9, 0.4, 0.2]
    return lander
end

@testset "Landing-target desirability" begin
    @testset "ranking utilities remain usable without branch utilities" begin
        standalone = joinpath(@__DIR__, "standalone_sort.jl")
        output = read(`$(Base.julia_cmd()) --startup-file=no $standalone`, String)
        @test strip(output) == "3,1,2"
    end

    @testset "inactive NaN metrics cannot change target ranking" begin
        lander = three_targets()
        sort_des_score!(lander)
        @test lander.λ_targs == [3, 1, 2]
    end

    @testset "inactive NaN metrics cannot suppress target switching" begin
        lander = three_targets()
        @test switch_decision(lander, 1) === true
    end

    @testset "inactive infinite metrics are ignored" begin
        lander = three_targets(; inactive=Inf)
        sort_des_score!(lander)
        @test lander.λ_targs == [3, 1, 2]
        @test switch_decision(lander, 1) === true
    end

    @testset "an active nonfinite metric is rejected" begin
        lander = three_targets()
        lander.p_targs["prox_clust"][1] = NaN
        @test_throws ArgumentError sort_des_score!(lander)
        @test_throws ArgumentError switch_decision(lander, 1)
    end

    @testset "changing active weights changes the relevant evidence" begin
        lander = three_targets()
        lander.w_des .= [1, 0, 0, 0, 0]
        lander.p_targs["pcd"] .= [0.2, 0.8, 0.4]
        lander.p_targs["prox_clust"] .= NaN  # Now inactive.
        @test switch_decision(lander, 2) === true
    end

    @testset "finite defaults preserve existing decisions" begin
        lander = three_targets(; inactive=1.0)
        sort_des_score!(lander)
        @test lander.λ_targs == [3, 1, 2]
        @test switch_decision(lander, 1) === true
        @test switch_decision(lander, 3) === false
    end
end
