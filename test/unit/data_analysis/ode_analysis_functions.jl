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

    # 3. Throw an error if the x_ss does not have the same length as the number of nodes
    # We test with a difficult case of a non-totally connected graph
    g = SimpleGraph(4)
    add_edge!(g, 1, 2)
    add_edge!(g, 3, 4)
    # Giant component size: 2
    # If x_ss has length 4, it should throw an error
    x_ss = [0.1, 0.2, 0.3, 0.4]
    @test_throws ArgumentError get_lambda(g, interaction_mode, taxation_mode, T, x_ss)

    # 4. Test get_lambda for a known case of a disconnected graph from our simulations
    # A test for solving ODE with disconnected networks
    # Graph parameters
    # Number of nodes
    N = 128
    # Mean degree
    k_mean = 8
    # Probability of connection
    p = k_mean/(N-1)
    seed_er = 1
    g = erdos_renyi(N,p,seed=seed_er)
    while is_connected(g)
        seed_er += 1
        g = erdos_renyi(N,p,seed=seed_er)
    end
    println("Seed for disconnected graph: ", seed_er)
    println("Number of connected components: ", length(connected_components(g)))
    N_gc = nv(get_giant_component(g))
    x_ss = solve_ode_net_SS(
        g,
        interaction_mode,
        taxation_mode,
        T,
        seed_er
        ).u
    @test size(get_lambda(g, interaction_mode, taxation_mode, T, x_ss)) == (N_gc, N_gc)
end
