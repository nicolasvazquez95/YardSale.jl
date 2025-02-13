using YardSale, Test, Graphs

# Test get_lambda and get_max_eigenvalue
@testset "Test get_lambda and get_max_eigenvalue" begin
    # 1. Parameters for the test (an example graph, actual simulation)
    N = 32
    k_mean = 8
    p = k_mean/(N-1)
    seed = 42
    g = erdos_renyi(N, p, seed=seed)
    T = 1.0
    interaction_mode = "A"
    taxation_mode = "A"
    # Run a simulation to get the steady state
    x_ss = solve_ode_net_SS(
        g,
        interaction_mode,
        taxation_mode,
        T,
        seed
        ).u
    # Get the lambda matrix
    Lambda = get_lambda(g, interaction_mode, taxation_mode, T, x_ss)
    # Get the maximum eigenvalue
    max_eigenvalue = get_max_eigenvalue(Lambda)
    # Test get_lambda
    @test size(Lambda) == (N,N)
    # Test get_max_eigenvalue
    @test typeof(max_eigenvalue) == Float64

    # 2. Test get_max_eigenvalue for a known matrix
    A = [0 1; 1 0]
    # Eigenvalues: 1, -1
    @test get_max_eigenvalue(A) ≈ 1.0
end
