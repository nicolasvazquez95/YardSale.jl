# Test functions for the Monte Carlo simulations.

using YardSale, Test

# parameters
N = 32
W_N = 1.0
# Uniform
w1 = mc_set_initial_conditions(N, W_N, "uniform")
# Random
w2 = mc_set_initial_conditions(N, W_N, "random")
w2_max = maximum(w2)
# Custom: totally condensed
w_custom = zeros(N)
w_custom[1] = 32.0
w3 = mc_set_initial_conditions(N, W_N, "custom", w_custom)

# Test get_x1
@testset "get_x1" begin
    @test get_x1(w1) ≈ 1/N
    @test get_x1(w2) ≈ w2_max/sum(w2)
    @test get_x1(w3) ≈ 1.0
end

@testset "get_xi" begin
    @test get_xi(w1, 1) == get_x1(w1)
    @test get_xi(w3, 1) == get_x1(w3)
    @test get_xi(w3, 2) ≈ 0.0
end

@testset "get_gini" begin
    @test get_gini(w1) ≈ 0.0
    @test 0.0 ≤ get_gini(w2) ≤ 1.0
    # The Gini coefficient of a totally condensed distribution is 1.0.
    # The formula is multiplied by (N / (N - 1)) to correct the bias.
    @test get_gini(w3) * (N / (N - 1)) ≈ 1.0
end

@testset "get_lorenz" begin
    @test get_lorenz(w1) ≈ cumsum(ones(N)) / N
    @test get_lorenz(w3) ≈ reverse(w3)/sum(w3)
end

@testset "get_R" begin
    @test get_R(w1) ≈ 0
    @test get_R(w3) ≈ 1
end

@testset "get_u" begin
    @test get_u(w1) ≈ 0
    @test get_u(w3) ≈ 1
end
