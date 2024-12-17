# ode_solvers.jl
# Here we define the functions to solve the ODEs for the wealth exchange models.
# We use the DifferentialEquations.jl package to solve the ODEs.
# To calculate the κ matrix and the β vector we use the functions defined in kappa_beta.jl.

# Derivatives
"""
    dxdt_net!(dxdt, x, p, t)
In-place function to calculate the derivatives of the network model.
# Arguments
    dxdt::Vector{Float64}: Vector of derivatives.
    x::Vector{Float64}: Vector of wealths.
    p::Tuple{Int64, Int64, Vector{Tuple{Int64, Int64}}, Vector{Float64}, Vector{Float64},
 Float64, Float64}: Tuple of parameters.
    t::Float64: Time.

# Details
The network model is given by the following set of ODEs:

```math
\\dot{x}_i = \\frac{T}{N} (- \\beta_i x_i + \\sum_j^N \beta_j x_j) + \\sum_j^N \\kappa_{ij}
(x_i - x_j) \\min(x_i, x_j)
```
where ``x_i`` is the wealth of node ``i``, ``T`` is temperature, ``N`` is
the number of nodes, ``beta_i`` is the redistribution parameter of node ``i``,
``kappa_{ij}`` is the exchange parameter between nodes ``i`` and ``j``.
"""
function dxdt_net!(dxdt,x,p,t)
    # Unpack parameters
    N, l, edgelist, kappa, beta, T_n, n_1 = p
    # Set the derivatives to zero
    fill!(dxdt, 0.0)
    # Calculate the derivatives
    @inbounds for link in 1:l
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

    @. dxdt += T_n * (bjxj_n - b*x)
end

# Solve ODE
"""
    solve_ode_net(g, tspan, integrator, IM, TM, T, seed; kwargs...)
Solve the ODE for the network model.
# Arguments
    g::SimpleGraph{<:Integer}: Graph.
    tspan::Tuple{Float64, Float64}: Tuple with initial and final time.
    integrator::Function: ODE solver.
    IM::String: Interaction matrix.
    TM::String: Topology matrix.
    T::Real: Temperature.
    seed::Integer: Random seed.
    kwargs...: Keyword arguments for the ODE solver.
# Details
This function solves the ODE for the network model. It calculates the kappa and beta
parameters, sets the initial conditions, and solves the ODE. It returns the solution.

# Example
```julia
using DifferentialEquations
using Graphs
using YardSale

# Create a graph
g = erdos_renyi_graph(100, 0.1, seed=42)
IM, TM = "A","A"
T = 1.0
seed = 42
tspan = (0.0, 10.0)
sol = solve_ode_net(g, tspan, Tsit5(), IM, TM, T, seed)
```
"""
function solve_ode_net(
    g::SimpleGraph{<:Integer},
    tspan::Tuple{Float64, Float64},
    integrator::Function,
    IM::String,
    TM::String,
    T::Real,
    seed::Integer;
    # Keyword arguments for the ODE solver
    kwargs...
    )

    # Set the random seed
    Random.seed!(seed)
    # Get the number of nodes
    N = nv(g)
    # Get the number of links
    l = length(edges(g))
    # Get the edgelist
    edgelist = collect(edges(g))

    # Calculate kappa and beta
    kappa, beta = get_kappa_beta(g, IM, TM)

    # Other constants
    n_1 = 1/N
    T_n = T/N

    # Parameters
    p = (N, l, edgelist, kappa, beta, T_n, n_1)

    # Initial conditions
    gauss = Normal(0.0, 0.01)
    # x_i = (1/N) * (1 + epsilon_i)
    x = n_1 * (ones(N) + rand(gauss, N))

    # Define the ODE problem
    prob = ODEProblem(dxdt_net!, x, tspan, p)

    # Solve the ODE
    sol = integrator(prob, integrator; kwargs...)

    # Check if the solver was successful
    if sol.retcode != :Success
        println("Solver failed with retcode: $(sol.retcode)")
    end
    return sol
end
