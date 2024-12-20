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
    gauss = Normal(0.0, W_N)
    noise = rand(gauss, N)
    custom_w = fill(W_N, N) + typeof(W_N).(noise)
    custom_w *= W/sum(custom_w)
    @test_nowarn begin
        w = mc_set_initial_conditions(N, W_N, "custom", custom_w)
    end
    @test sum(w) ≈ W
end
