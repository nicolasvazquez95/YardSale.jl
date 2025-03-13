# ode_solvers.jl
# Here we define the functions to solve the ODEs for the wealth exchange models.
# We use the DifferentialEquations.jl package to solve the ODEs.
# To calculate the κ matrix and the β vector we use the functions defined in kappa_beta.jl.

# Derivatives
"""
    dxdt_net!(dxdt, x, p, t)
In-place function to calculate the derivatives of the network model.
# Arguments
    dxdt::Vector{<:Real}: Vector of derivatives.
    x::Vector{<:Real}: Vector of wealths.
    p::Tuple{<:Integer, <:Integer, Vector{Vector{<:Integer}},
    Vector{<:Integer}, Vector{<:Array}, ::Real, ::Real}: Parameters.
    t::Float64: Time.

# Details
The network model is given by the following set of ODEs:

```math
\\dot{x}_i = \\frac{T}{N} (- \\beta_i x_i + \\frac{1}{N}\\sum_j^N \\beta_j x_j) + \\sum_j^N \\kappa_{ij}
(x_i - x_j) \\min(x_i, x_j)
```
where ``x_i`` is the wealth of node ``i``, ``T`` is temperature, ``N`` is
the number of nodes, ``\\beta_i`` is the redistribution parameter of node ``i``,
``\\kappa_{ij}`` is the exchange parameter between nodes ``i`` and ``j``.
"""
function dxdt_net!(dxdt,x,p,t)
    # Unpack parameters
    N, l, edgelist, kappa, beta, T_n, n_1 = p
    # Set the derivatives to zero
    fill!(dxdt, 0.0)
    # Calculate the derivatives
    for link in 1:l
        # Get nodes, wealths, and kappa
        i, j = edgelist[link]
        xi, xj = x[i], x[j]
        kappa_ij = kappa[link]
        # Calculate the exchange term
        exch_ij = kappa_ij * (xi - xj) * min(xi, xj)
        # Update the derivatives
        dxdt[i] += exch_ij
        dxdt[j] -= exch_ij
    end
    # Add the redistribution term
    bjxj_n = n_1 * sum(beta .* x)

    @. dxdt += T_n * (bjxj_n - beta*x)
end

# Solve ODE
"""
    solve_ode_net(g, tspan, integrator, interaction_mode, taxation_mode, T, seed; kwargs...)
Solves the ODE for the network model.
# Arguments
    g::SimpleGraph{<:Integer}: Graph.
    tspan::Tuple{Float64, Float64}: Tuple with initial and final time.
    interaction_mode::String: Interaction mode.
    taxation_mode::String: Taxation mode.
    T::Real: Temperature.
    seed::Integer: Random seed.
# Optional arguments
    `integrator::SciMLAlgorithm`: An instante of the integrator to use. Default is Tsit5().
    initial_conditions::String=nothing: Initial condition. Options are
    "uniform", "random", "noisy" and "custom". Default is "uniform". If "custom" is
    chosen, the x0 argument must be provided.
    x0::Union{Nothing, Vector{<:Real}}=nothing: Initial wealth distribution.
    `kwargs...`: Additional arguments for the solver.
# Details
This function solves the ODE for the network model. It calculates the kappa and beta
parameters, sets the initial conditions, and solves the ODE. It returns the solution in
the standard `DifferentialEquations.jl` format.
A brief description of the ODE approximation model can be found in the `dxdt_net!` function.
The function uses the `DifferentialEquations.jl` package to solve the ODE.
All the parameters for the solver can be passed as keyword arguments.

The initial conditions can be set to "noisy", "random", "uniform", or a custom vector.
- `"noisy"`: Initial conditions are set to a random value around 1/N.
``x = (1/N) * (1 + \\epsilon)`` where `ϵ` is white noise with μ=0 and σ=0.01.
- `"random"`: Initial conditions are set to a random value, normalized to sum to 1.
- `"uniform"`: Initial conditions are set to 1/N.
- `"custom"`: Initial conditions are set to a custom vector. The vector must sum to 1.
Also, ``x_i`` must be positive for all ``i``.

# Returns
    sol::ODESolution: Solution of the ODE.
# Examples
```julia
using DifferentialEquations, Graphs, YardSale
# Create a graph
g = erdos_renyi(100, 0.1, seed=42)
interaction_mode, taxation_mode = "A","A"
T = 1.0
seed = 42
tspan = (0.0, 10.0)
sol1 = solve_ode_net(g, tspan, interaction_mode, taxation_mode, T, seed)
sol2 = solve_ode_net(g, tspan, interaction_mode, taxation_mode, T, seed; integrator=RK4())
sol3 = solve_ode_net(g, tspan, interaction_mode, taxation_mode, T, seed;
integrator=RK4(), reltol=1e-6, abstol=1e-6)
```
"""
function solve_ode_net(
    g::SimpleGraph{<:Integer},
    tspan::Tuple{<:Real, <:Real},
    interaction_mode::String,
    taxation_mode::String,
    T::Real,
    seed::Integer;
    integrator::SciMLAlgorithm = Tsit5(),
    initial_conditions::Union{String, Vector{<:Real}}="noisy",
    x0::Union{Nothing, Vector{<:Real}}=nothing,
    kwargs...
    )

    # Set the random seed
    Random.seed!(seed)
    # Check if the graph is connected
    if !is_connected(g)
        @warn "The graph is not connected. Using the giant component for the simulation."
    end
    # Get the graph data
    N, l , edgelist, gc = get_graph_data(g)
    # Calculate kappa and beta
    kappa, beta = get_kappa_beta(gc, interaction_mode, taxation_mode)
    # Other constants
    n_1 = 1/N
    T_n = T/N
    # Parameters
    p = (N, l, edgelist, kappa, beta, T_n, n_1)
    # Set initial conditions
    x = ode_set_initial_conditions(N, initial_conditions, x0)

    # Define the ODE problem
    prob = ODEProblem(dxdt_net!, x, tspan, p)

    # Solve the ODE
    sol = solve(prob, integrator; kwargs...)

    # Check if the solver was successful
    if sol.retcode != :Success
        println("Solver failed with retcode: $(sol.retcode)")
    end
    return sol
end

"""
    solve_ode_net_SS(g, interaction_mode, taxation_mode, T, seed;
    initial_conditions="noisy",integrator=Tsit5(), kwargs...)
Solves the ODE for the network model using a steady state solver.
# Arguments
    g::SimpleGraph{<:Integer}: Graph.
    interaction_mode::String: Interaction mode.
    taxation_mode::String: Taxation mode.
    T::Real: Temperature.
    seed::Integer: Random seed.

# Optional arguments
    `initial_conditions::Union{String, Vector{<:Real}}`="noisy": Initial conditions
    for the ODE. Default is "noisy", which sets the initial conditions to a random value
    around 1/N. Other options are "random" and "uniform". If a vector is passed, it will be
    used as the initial conditions.
    `x0::Union{Nothing, Vector{<:Real}}`=nothing: Initial wealth distribution.
    `integrator::SciMLAlgorithm`: An instance of the integrator to use. Default is Tsit5().
    `kwargs...`: Additional arguments for the solver.
# Details
This function solves the ODE for the network model using a steady state solver.
It is similar to `solve_ode_net`, but it uses the `DynamicSS` solver to find
the steady state.
Instead of an ODE problem, it uses a `SteadyStateProblem` to solve the ODE.
Reference: https://docs.sciml.ai/DiffEqDocs/stable/types/steady_state_types/

As in `solve_ode_net`, the initial conditions can be set to "noisy", "random", "uniform",
or a custom vector.
# Returns
    sol::ODESolution: Solution of the ODE.
# Examples
```julia
using DifferentialEquations, Graphs, YardSale
# Create a graph
g = erdos_renyi(100, 0.1, seed=42)
interaction_mode, taxation_mode = "A","A"
T = 1.0
seed = 42
sol1 = solve_ode_net_SS(g, interaction_mode, taxation_mode, T, seed)
sol2 = solve_ode_net_SS(g, interaction_mode, taxation_mode, T, seed; integrator=RK4())
sol3 = solve_ode_net_SS(g, interaction_mode, taxation_mode, T, seed;
                            integrator=RK4(),reltol=1e-6, abstol=1e-6
                            )
```
"""
function solve_ode_net_SS(
    g::SimpleGraph{<:Integer},
    interaction_mode::String,
    taxation_mode::String,
    T::Real,
    seed::Integer;
    initial_conditions::Union{String, Vector{<:Real}}="noisy",
    x0::Union{Nothing, Vector{<:Real}}=nothing,
    integrator::SciMLAlgorithm = Tsit5(),
    kwargs...
    )
    # Set the random seed
    Random.seed!(seed)
    # Check if the graph is connected
    if !is_connected(g)
        @warn "The graph is not connected. Using the giant component for the simulation."
    end
    # Get the graph data
    N, l , edgelist, gc = get_graph_data(g)
    # Calculate kappa and beta
    kappa, beta = get_kappa_beta(gc, interaction_mode, taxation_mode)

    # Other constants
    n_1 = 1/N
    T_n = T/N

    # Parameters
    p = (N, l, edgelist, kappa, beta, T_n, n_1)

    # Set initial conditions
    x = ode_set_initial_conditions(N, initial_conditions, x0)

    # Define the ODE problem

    # Set the time span to an arbitrary large value, only to define the initial ODE problem
    tspan = (0.0, 1e6)
    prob = ODEProblem(dxdt_net!, x, tspan, p)

    # Define the Steady State problem from the ODE problem
    ss_prob = SteadyStateProblem(prob)

    # Solve the ODE with a steady state solver (DynamicSS)
    sol = solve(ss_prob, DynamicSS(integrator); kwargs...)

    # Check if the solver was successful
    if !SciMLBase.successful_retcode(sol)
        println("Solver failed with retcode: $(sol.retcode)")
        return
    end
    return sol
end
