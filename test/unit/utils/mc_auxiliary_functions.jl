# Test mc_auxiliary_functions

using YardSale, Test, Distributions

# Define main parameters
N = 32
W_N = 1.0
W = W_N * N

# Test mc_set_initial_conditions
@testset "mc_set_initial_conditions" begin
    # Test uniform initial conditions
    w = mc_set_initial_conditions(N, W_N, "uniform")
    @test w ≈ fill(W_N, N)

    # Test random initial conditions
    w = mc_set_initial_conditions(N, W_N, "random")
    @test sum(w) ≈ W

    # Test noisy initial conditions
    w = mc_set_initial_conditions(N, W_N, "noisy")
    @test sum(w) ≈ W

    # Test custom initial conditions
    custom_w = fill(0.5, N)
    # 1. It should throw a warning
    @test_logs (:warn,"The sum of w is different from W. It will be normalized to W by default.") begin
        w = mc_set_initial_conditions(N, W_N, "custom", custom_w)
    end
    @test sum(w) ≈ W

    # 2. It should not throw a warning
    custom_w = rand(typeof(W_N),N)
    custom_w *= W/sum(custom_w)
    @test_nowarn begin
        w = mc_set_initial_conditions(N, W_N, "custom", custom_w)
    end
    @test sum(w) ≈ W

    # 3. Errors
    # 3.1. Invalid initial condition
    @test_throws ArgumentError mc_set_initial_conditions(N, W_N, "invalid")
    # 3.2. Custom initial conditions without w
    @test_throws ArgumentError mc_set_initial_conditions(N, W_N, "custom")
    # 3.3. Custom initial conditions with wrong length of w
    @test_throws ArgumentError mc_set_initial_conditions(N, W_N, "custom", fill(0.5, N+1))
    # 3.4. Custom initial conditions with negative elements
    @test_throws ArgumentError mc_set_initial_conditions(N, W_N, "custom", fill(-0.5, N))
end

@testset "mc_set_beta" begin
    # Test beta0
    beta = mc_set_beta(N, "uniform", 1.0)
    @test beta == fill(1.0, N)

    beta = mc_set_beta(N, "random")
    @test length(beta) == N

    beta = mc_set_beta(N, "noisy", 0.8)
    @test length(beta) == N

    beta0 = rand(typeof(W_N),N)
    beta = mc_set_beta(N, "custom", beta0)
    @test beta0 == beta
end
