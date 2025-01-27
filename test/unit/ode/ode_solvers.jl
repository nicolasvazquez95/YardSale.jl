# Test file for the ODE solvers in the ode_solvers.jl file.
using YardSale, Test, Graphs

@testset "ODE in disconnected networks" begin
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
    println("Number of connected components: ", connected_components(g))

    # Size of giant component
    n_gc = nv(get_giant_component(g))

    # Interaction/Taxation Modes
    interaction_mode = "A"
    taxation_mode = "B"
    # Initial conditions
    initial_conditions = "noisy"
    T = 1.0
    seed = 42
    sol = solve_ode_net_SS(g, interaction_mode, taxation_mode, T, seed, initial_conditions)

    # Test that the solution has the correct size
    @test length(sol) == n_gc
end
