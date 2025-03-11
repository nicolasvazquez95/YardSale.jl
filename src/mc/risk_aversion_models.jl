"""
    YS_base_risk(N, W_N, T, initial_conditions, beta; seed = nothing)

Run the Yard-Sale model with risk aversion for N agents with wealth exchange
in a fully connected network.

# Arguments
    N::Integer: Number of agents
    W_N::Real: Wealth/N quotient
    steps::Integer: Number of MC steps
    initial_conditions::AbstractVector{<:Real}: Initial conditions
    beta::AbstractVector{<:Real}: Beta values (risk aversion)
    w0::Union{Nothing, Vector{<:Real}}=nothing: Initial conditions (optional)

    seed::Union{Nothing, Int64}: Seed for reproducibility (optional)

# Details
The model is defined by the following wealth exchange rule:
```math
\\Delta w_{ij} = \\min((1-\\beta_i) w_i, (1 - \\beta_j) w_j)
```
where ``\\beta_i`` is the risk aversion of agent i.

# Output
    w_t::Matrix{Float32}: T x N matrix with the wealth of the agents at each time step
"""
function YS_base_risk(
    N::Integer, # Number of agents
    W_N::Real, # Wealth/N quotient
    steps::Integer, # Number of MC steps
    beta::String, # Beta values (risk aversion)
    seed::Integer; # Seed for reproducibility
    initial_conditions::String="uniform",
    w0::Union{Nothing, Vector{<:Real}}=nothing, # Initial conditions (optional)
    beta0::Union{Nothing, Real, Vector{<:Real}}=0.95, # Initial beta values (optional)
    save_every::Union{Nothing, Integer}=nothing
    )

    # Initialize the model with the initial conditions
    w = mc_set_initial_conditions(N, W_N, initial_conditions, w0)
    W = sum(w) # Total wealth

    # Initialize the beta values
    beta = mc_set_beta(N, beta, beta0)
    # Define risk propension as r = 1 - beta
    r = 1 .- beta # Risk propension (it is computationally more efficient to calculate it once)

    # Initialize the wealth matrix to store the results
    ## Check the save_every argument
    if isnothing(save_every)
        save_every = N
    end
    ## Initialize the time series
    ### Each row of w_t is a checkpoint
    w_t = zeros(typeof(W_N), 1 +(steps ÷ save_every), N)
    ## Save the initial condition
    w_t[1, :] .= w
    # Index to iterate over w_t
    idx = 2

    # Initialize the random number generator
    if isnothing(seed) == false
        Random.seed!(seed)
    end

    # Define some useful constants
    n_2 = N ÷ 2 # For the MC step definition

    # Run the model
    for t in 2:steps
        # Random numbers for the mc step
        random_numbers = rand(typeof(W_N), n_2) # n_2 random numbers for the MC step. It is more efficient to calculate them once.
        eta = 2 * (random_numbers .> 0.5) .-1 # eta = 1 if random_number > 0.5, eta = -1 otherwise
        for exch in 1:n_2
            # Choose two agents at random
            i, j = rand(1:N, 2)
            # If i==j, choose another agent
            while i == j
                i, j = rand(1:N, 2)
            end
            # Calculate the wealth exchange
            w_i, w_j = w[i], w[j] # Get the wealth of the agents
            r_i, r_j = r[i], r[j] # Get the risk propension of the agents

            dw_exch = eta[exch] * Δw_risk(w_i, w_j, r_i, r_j) # Calculate the wealth exchange
            # Update the wealth of the agents
            w[i] += dw_exch
            w[j] -= dw_exch
        end
        # Store the wealth of the agents
        if t % save_every == 0
            # Re normalize the wealth
            w ./= sum(w)
            w .*= W

            w_t[idx, :] .= w
            idx += 1
            # Check for negative wealth
            if any(w .< 0)
                throw(ArgumentError("Negative wealth detected. Simulation stopped."))
                return w_t
            end
        end
    end
    return w_t
end

"""
    YS_net_risk(g, W_N, T, initial_conditions, beta; seed = nothing, exchange_mode)

Run the Yard-Sale model with risk aversion for N agents with wealth exchange
in a complex network.

# Arguments
    g::SimpleGraph{<:Integer}: Graph
    W_N::Real: Wealth/N quotient
    steps::Integer: Number of MC steps
    initial_conditions::AbstractVector{<:Real}: Initial conditions
    beta::AbstractVector{<:Real}: Beta values (risk aversion)
    w0::Union{Nothing, Vector{<:Real}}=nothing: Initial conditions (optional)
    seed::Union{Nothing, Int64}: Seed for reproducibility (optional)
    exchange_mode::String: "link" or "node"

# Details
The model is defined by the following wealth exchange rule:
```math
\\Delta w_{ij} = \\min((1-\\beta_i) w_i, (1 - \\beta_j) w_j)
```
where ``\\beta_i`` is the risk aversion of agent i.

At each time step, two agents are chosen at random. The wealth exchange is calculated
according to the rule above, and the wealth of the agents is updated.

Two exchange modes are available:
- "link": One link chosen at random, and the wealth exchange is calculated between the two
agents.
- "node": A node is chosen at random, and one of its neighbors is chosen at random.

# Output
    w_t::Matrix{Float32}: T x N matrix with the wealth of the agents at each time step
"""
function YS_net_risk(
    g::SimpleGraph{<:Integer}, # Graph
    W_N::Real, # Wealth/N quotient
    steps::Integer, # Number of MC steps
    beta::String, # Beta values (risk aversion)
    seed::Integer; # Seed for reproducibility
    initial_conditions::String="uniform",
    w0::Union{Nothing, Vector{<:Real}}=nothing, # Initial conditions (optional)
    beta0::Union{Nothing, Real, Vector{<:Real}}=0.95, # Initial beta values (optional)
    save_every::Union{Nothing, Integer}=nothing
    )

    N = nv(g) # Number of agents
    edgelist = [(e.src,e.dst) for e in edges(g)]
    neighbors_list = [neighbors(g, i) for i in 1:N]
    if exchange_mode ∉ ["link", "node"]
        error("exchange_mode must be 'link' or 'node'")
    end

    # Initialize the model with the initial conditions
    w = mc_set_initial_conditions(N, W_N, initial_conditions, w0)
    W = sum(w) # Total wealth

    # Initialize the beta values
    beta = mc_set_beta(N, beta, beta0)

    # Define risk propension as r = 1 - beta
    r = 1 .- beta # Risk propension (it is computationally more efficient to calculate it once)

    # Initialize the wealth matrix to store the results
    ## Check the save_every argument
    if isnothing(save_every)
        save_every = N
    end
    ## Initialize the time series
    ### Each row of w_t is a checkpoint
    w_t = zeros(typeof(W_N), 1 +(steps ÷ save_every), N)
    ## Save the initial condition
    w_t[1, :] .= w
    # Index to iterate over w_t
    idx = 2

    # Initialize the random number generator
    if isnothing(seed) == false
        Random.seed!(seed)
    end

    # Define some useful constants
    n_2 = N ÷ 2 # For the MC step definition

    # Run the model
    for t in 2:steps
        # Random numbers for the mc step
        random_numbers = rand(typeof(W_N), n_2) # n_2 random numbers for the MC step. It is more efficient to calculate them once.
        eta = 2 * (random_numbers .> 0.5) .-1 # eta = 1 if random_number > 0.5, eta = -1 otherwise
        if exchange_mode == "link"
            links = rand(edgelist, n_2)
            for exch in 1:n_2
                # Choose two agents at random
                i, j = links[exch]

                # Calculate the wealth exchange
                w_i, w_j = w[i], w[j] # Get the wealth of the agents
                r_i, r_j = r[i], r[j] # Get the risk propension of the agents

                dw_exch = eta[exch] * Δw_risk(w_i, w_j, r_i, r_j) # Calculate the wealth exchange
                # Update the wealth of the agents
                w[i] += dw_exch
                w[j] -= dw_exch
            end
            # Store the wealth of the agents
            if t % save_every == 0
                # Re normalize the wealth
                w ./= sum(w)
                w .*= W

                w_t[idx, :] .= w
                idx += 1
                # Check for negative wealth
                if any(w .< 0)
                    throw(ArgumentError("Negative wealth detected. Simulation stopped."))
                    return w_t
                end
            end
        elseif exchange_mode == "node"
            nodes_i = rand(1:N, n_2)
            nodes_j = [rand(neighbors_list[i]) for i in nodes_i]
            for exch in 1:n_2
                # Choose two agents at random
                i, j = nodes_i[exch], nodes_j[exch]

                # Calculate the wealth exchange
                w_i, w_j = w[i], w[j] # Get the wealth of the agents
                r_i, r_j = r[i], r[j] # Get the risk propension of the agents

                dw_exch = eta[exch] * dw(w_i, w_j, r_i, r_j) # Calculate the wealth exchange
                # Update the wealth of the agents
                w[i] += dw_exch
                w[j] -= dw_exch
            end
            # Store the wealth of the agents
            if t % save_every == 0
                # Re normalize the wealth
                w ./= sum(w)
                w .*= W

                w_t[idx, :] .= w
                idx += 1
                # Check for negative wealth
                if any(w .< 0)
                    throw(ArgumentError("Negative wealth detected. Simulation stopped."))
                    return w_t
                end
            end
        end
    end
    return w_t
end
