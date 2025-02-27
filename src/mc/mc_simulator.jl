"""
    EYSM_base_full(N, W_N, chi, zeta, f, steps,seed; w=nothing,initial_conditions="uniform",
    save_every=nothing)
Runs a Monte Carlo simulation of the Extended Yard-Sale Model and returns the
whole time series of the wealth distribution.
# Arguments
    N::Integer: Number of agents.
    W_N::Real: Mean wealth of the agents. It is the total wealth divided by the number
    of agents.
    chi::Real: Taxation rate.
    zeta::Real: Wealth-Attained-Advantage parameter.
    f::Real: Fraction of the wealth that is redistributed.
    steps::Integer: Number of steps of the simulation, measured in Monte Carlo steps.
    seed::Integer: Seed for the random number generator.
# Optional arguments
    w0::Union{Nothing, Vector{<:Real}}=nothing: Initial wealth distribution.
    Only used if initial_conditions="custom". Default is nothing.
    initial_conditions::String=nothing: Initial condition. Options are
    "uniform", "random", "noisy" and "custom". Default is "uniform". If "custom" is
    chosen, the w0 argument must be provided.
    save_every::Union{Nothing, Integer}=nothing: Save the wealth distribution every
    save_every steps. Default is nothing, which means saving every N steps.
# Returns

# Details
This function tries to reproduce as closely as possible the results of the original EYSM
proposed by Boghosian et al. in 2017. The model has two parameters, and it is defined for
a fully connected network of ``N`` agentes, each one with a wealth ``w_i``, such that the
total wealth ``\\sum_i w_i = W`` is conserved.

The interaction rule is defined as follows:
1. Randomly select two agents i and j.
2. The amount of wealth exchanged is given by a fraction of the minimum between
the wealths of the two agents.
```math
\\Delta w = f\\min(w_i, w_j)\\eta_{ij}
```
where ``\\eta_{ij}`` is a stochastic variable with values -1 or 1. The expected value of
``\\eta_{ij}`` is biased towards the richest agent, and it is given by
```math
\\langle\\eta_{ij}\\rangle = \\zeta \\frac{w_i - w_j}{W}
```
3. The wealth of the agents is updated as follows:
```math
w_i \\to w_i - \\Delta w + \\chi (\\frac{W}{N} - w_i)
w_j \\to w_j + \\Delta w + \\chi (\\frac{W}{N} - w_j)
````
where ``\\chi`` represents the taxation and redistribution rate.

In the original paper, the authors work analitically, using a Fokker-Planck approach. In this
implementation, we use a Monte Carlo simulation to obtain similar results.

# Examples
```julia
w_t_1 = EYSM_base_full(1000, 1, 0.1, 0.1, 1, 1000, 42)
w_t_2 = EYSM_base_full(64, 1.0f0, 0.1f0, 0.1f0, 1, 1000, 42, initial_conditions="random")
w_t_3 = EYSM_base_full(32, 1, 0.1, 0.1, 1, 1000, 42, initial_conditions="custom", w=rand(1000))
```
"""
function EYSM_base_full(
    N::Integer,
    W_N::Real,
    chi::Real,
    zeta::Real,
    f::Real,
    steps::Integer,
    seed::Integer;
    w0::Union{Nothing, Vector{<:Real}}=nothing,
    initial_conditions::String="uniform",
    save_every::Union{Nothing, Integer}=nothing
)
    # Set the seed
    Random.seed!(seed)
    # Initialize wealth distribution
    w = mc_set_initial_conditions(N, W_N, initial_conditions, w0)
    # Calculate the total wealth
    W = sum(w)
    # Calculate some useful constants
    chif_N = chi * f / N
    zeta_W = zeta / W

    # Initialize the time series
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

    # Simulation loop
    for t in 1:steps
        for exch in 1:N
            # Select two agents
            i, j = rand(1:N, 2)
            # Avoid self-exchanges
            while i == j
                j = rand(1:N)
            end
            # Calculate the wealth exchange
            wi, wj = w[i], w[j]
            δw = Δw(f, wi, wj, zeta_W)
            # Calculate the redistribution terms
            redist_i = EYSM_base_redistribution(wi, W_N, chif_N)
            redist_j = EYSM_base_redistribution(wj, W_N, chif_N)
            # Update the wealth
            w[i] = wi + δw + redist_i
            w[j] = wj - δw + redist_j
        end
        # Save the wealth distribution
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
    EYSM_net_full(g::SimpleGraph,W_N, interaction_mode, taxation_mode, chi, zeta, f, steps,
    seed; w=nothing,initial_conditions="uniform",save_every=nothing)
Runs a Monte Carlo simulation of the Extender Yard-Sale model on an arbitrary network and
returns the whole time series of the wealth distribution.
# Arguments
    g::SimpleGraph{<:Integer}: Graph representing the network of agents.
    W_N::Real: Mean wealth of the agents. It is the total wealth divided by the number
    of agents.
    interaction_mode::String: Interaction mode. Options are "A" and "B".
    taxation_mode::String: Taxation mode. Options are "A" and "B".
    chi::Real: Taxation rate.
    zeta::Real: Wealth-Attained-Advantage parameter.
    f::Real: Fraction of the wealth that is redistributed.
    steps::Integer: Number of steps of the simulation, measured in Monte Carlo steps.
    seed::Integer: Seed for the random number generator.
# Optional arguments
    w0::Union{Nothing, Vector{<:Real}}=nothing: Initial wealth distribution.
    Only used if initial_conditions="custom". Default is nothing.
    initial_conditions::String=nothing: Initial condition. Options are
    "uniform", "random", "noisy" and "custom". Default is "uniform". If "custom" is
    chosen, the w argument must be provided.
    save_every::Union{Nothing, Integer}=nothing: Save the wealth distribution every
    save_every steps. Default is nothing, which means saving every N steps.
# Details
This function is the generalization of the EYSM_base_full function to arbitrary networks.
In the extension, an ambiguity arises in the definition of the interaction rule and the
taxation rule. That is why we have to define previously the way in which agents will be selected
to exchange wealth and the way in which the taxation will be applied. For each of these two
steps, we have two options, which are defined by the interaction_mode and taxation_mode.
The interaction_mode defines the way in which agents will be selected to exchange wealth.
The options are:
- "A": Randomly select a link.
- "B": Randomly select a node and one of its neighbors.
The taxation_mode defines the way in which the taxation will be applied. The options are:
- "A": Tax exchanging agents.
- "B": Tax two random agents.
The equation governing the wealth exchange is a generalization of the one used in the
EYSM_base_full function:
```math
\\Delta w_i(t) = \\chi\\Delta t(-b_i(t) w_i(t) + \\frac{1}{N}\\sum_j b_j(t)w_j(t))
+ \\sum_j c_{ij}(t)\\sqrt{\\gamma\\Delta t} \\eta_{ij}(t)\\min(w_i(t), w_j(t))
```
where now ``b_i(t)`` represents an event of taxation over agent ``i`` and ``c_{ij}(t)``
represents an event of wealth exchange between agents ``i`` and ``j``.
To see the ODE system that describes an equivalent dynamics, see the documentation of the
`solve_ode_net` function.
# Returns
    w_t::Matrix{Real}: Time series of the wealth distribution. Each row is a checkpoint.
"""
function EYSM_net_full(
    g::SimpleGraph{<:Integer},
    W_N::Real,
    interaction_mode::String,
    taxation_mode::String,
    chi::Real,
    zeta::Real,
    f::Real,
    steps::Integer,
    seed::Integer;
    w0::Union{Nothing, Vector{<:Real}}=nothing,
    initial_conditions::String="uniform",
    save_every::Union{Nothing, Integer}=nothing
)
    # Get data from the graph
    N, l, edgelist, gc = get_graph_data(g)
    # Set the seed
    Random.seed!(seed)
    # Initialize wealth distribution
    w = mc_set_initial_conditions(N, W_N, initial_conditions, w0)
    # Calculate the total wealth
    W = sum(w)
    # Calculate some useful constants
    chif_N = chi * f / N
    chif_N2 = chi * f / (N^2)
    zeta_W = zeta / W

    # Initialize the time series
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

    # Best ways to access the graph:
    # Random nodes: Instantiate 1:N
    nodes = 1:N
    # Neighbors of node: Neighbors list or neighbors(g, node) from Graphs.jl
    nl = [neighbors(gc, i) for i in nodes]

    # Check
    if interaction_mode ∉ ["A","B"]
        throw(ArgumentError("Invalid interaction mode. Choose 'A' or 'B'."))
    elseif taxation_mode ∉ ["A","B"]
        throw(ArgumentError("Invalid taxation mode. Choose 'A' or 'B'."))
    end

    # Simulation loop
    # IMA: Randomly select a link
    if interaction_mode == "A"
        # TMA: Tax exchanging agents
        if taxation_mode == "A"
            for t in 1:steps
                for exch in nodes
                    i,j = rand(edgelist)
                    wi, wj = w[i], w[j]
                    # Calculate the wealth exchange
                    δw = Δw(f, wi, wj, zeta_W)
                    # Update the wealth and apply the taxation
                    w[i] += δw - chif_N * wi
                    w[j] += -δw - chif_N * wj
                    # Everybody receives the benefits of the taxation
                    # Broadcasting in TMA is more efficient! (I measured it)
                    @. w += chif_N2 * (wi + wj)
                end
                # Save the wealth distribution
                if t % save_every == 0
                    # Re normalize the wealth
                    w ./= sum(w)
                    w .*= W

                    @. w_t[idx, :] = w
                    idx += 1
                    # Check for negative wealth
                    if any(w .< 0)
                        throw(ArgumentError("Negative wealth detected. Simulation stopped."))
                        return w_t
                    end
                end
            end
        # TMB: Tax two random agents
        elseif taxation_mode == "B"
            for t in 1:steps
                taxed_agents = [(rand(nodes), rand(nodes)) for _ in 1:N]
                for exch in nodes
                    # Exchange agents
                    i,j = rand(edgelist)
                    wi, wj = w[i], w[j]
                    # Taxed agents
                    i_taxed, j_taxed = taxed_agents[exch]
                    wi_taxed, wj_taxed = w[i_taxed], w[j_taxed]
                    # Calculate the wealth exchange
                    δw = Δw(f, wi, wj, zeta_W)
                    # Update the wealths
                    w[i] += δw
                    w[j] -= δw
                    # Tax the agents
                    w[i_taxed] -= chif_N * wi_taxed
                    w[j_taxed] -= chif_N * wj_taxed
                    # Everybody receives the benefits of the taxation
                    # Trying to reduce the number of memory allocations by using a for loop
                    # It has a little impact on performance on TMB (I measured it)
                    tax_benefit = chif_N2 * (wi_taxed + wj_taxed)
                    for i in nodes
                        w[i] += tax_benefit
                    end
                end
                # Save the wealth distribution
                if t % save_every == 0
                    # Re normalize the wealth
                    w ./= sum(w)
                    w .*= W

                    @. w_t[idx, :] = w
                    idx += 1
                    # Check for negative wealth
                    if any(w .< 0)
                        throw(ArgumentError("Negative wealth detected. Simulation stopped."))
                        return w_t
                    end
                end
            end
        end
    # IMB: Randomly select a node and one of its neighbors
    elseif interaction_mode == "B"
        # TMA: Tax exchanging agents
        if taxation_mode == "A"
            for t in 1:steps
                rand_nodes = rand(1:N, N)
                rand_neighbors = [rand(nl[i]) for i in rand_nodes]
                for exch in nodes
                    i = rand_nodes[exch]
                    j = rand_neighbors[exch]
                    wi, wj = w[i], w[j]
                    # Calculate the wealth exchange
                    δw = Δw(f, wi, wj, zeta_W)
                    # Update the wealth
                    w[i] += δw - chif_N * wi
                    w[j] += -δw - chif_N * wj
                    # Everybody receives the benefits of the taxation
                    # Broadcasting in TMA is more efficient! (I measured it)
                    @. w += chif_N2 * (wi + wj)
                end
                # Save the wealth distribution
                if t % save_every == 0
                    # Re normalize the wealth
                    w ./= sum(w)
                    w .*= W

                    @. w_t[idx, :] = w
                    idx += 1
                    # Check for negative wealth
                    if any(w .< 0)
                        throw(ArgumentError("Negative wealth detected. Simulation stopped."))
                        return w_t
                    end
                end
            end
        # TMB: Tax two random agents
        elseif taxation_mode == "B"
            for t in 1:steps
                taxed_agents = [(rand(nodes), rand(nodes)) for _ in 1:N]
                rand_nodes = rand(1:N, N)
                rand_neighbors = [rand(nl[i]) for i in rand_nodes]
                for exch in nodes
                    # Exchange agents
                    i = rand_nodes[exch]
                    j = rand_neighbors[exch]
                    wi, wj = w[i], w[j]
                    # Taxed agents
                    i_taxed, j_taxed = taxed_agents[exch]
                    wi_taxed, wj_taxed = w[i_taxed], w[j_taxed]
                    # Calculate the wealth exchange
                    δw = Δw(f, wi, wj, zeta_W)
                    # Update the wealth
                    w[i] += δw
                    w[j] -= δw
                    # Tax the agents
                    w[i_taxed] -= chif_N * wi_taxed
                    w[j_taxed] -= chif_N * wj_taxed
                    # Everybody receives the benefits of the taxation
                    # Trying to reduce the number of memory allocations by using a for loop
                    # It has a little impact on performance on TMB (I measured it)
                    tax_benefit = chif_N2 * (wi_taxed + wj_taxed)
                    for i in nodes
                        w[i] += tax_benefit
                    end
                end
                # Save the wealth distribution
                if t % save_every == 0
                    # Re normalize the wealth
                    w ./= sum(w)
                    w .*= W

                    @. w_t[idx, :] = w
                    idx += 1
                    # Check for negative wealth
                    if any(w .< 0)
                        throw(ArgumentError("Negative wealth detected. Simulation stopped."))
                        return w_t
                    end
                end
            end
        end
    end
    return w_t
end

"""
    EYSM_base_callbacks(N, W_N, chi, zeta, f, steps, seed; w=nothing,
    initial_conditions="uniform",save_every=nothing, callbacks=nothing,
    callbacks_only=false)
Runs a Monte Carlo simulation of the Extended Yard-Sale Model on a fully connected network
and returns the whole time series of the wealth distribution and
the results of the callbacks.

# Arguments
    N::Integer: Number of agents.
    W_N::Real: Mean wealth of the agents. It is the total wealth divided by the number
    of agents.
    chi::Real: Taxation rate.
    zeta::Real: Wealth-Attained-Advantage parameter.
    f::Real: Fraction of the wealth that is redistributed.
    steps::Integer: Number of steps of the simulation, measured in Monte Carlo steps.
    seed::Integer: Seed for the random number generator.
# Optional arguments
    w0::Union{Nothing, Vector{<:Real}}=nothing: Initial wealth distribution.
    Only used if initial_conditions="custom". Default is nothing.
    initial_conditions::String=nothing: Initial condition. Options are
    "uniform", "random", "noisy" and "custom". Default is "uniform". If "custom" is
    chosen, the w argument must be provided.
    save_every::Union{Nothing, Integer}=nothing: Save the wealth distribution every
    save_every steps. Default is nothing, which means saving every N steps.
    callbacks::Union{Nothing, Dict{Symbol, Function}}=nothing: Dictionary of callbacks.
    The keys are the names of the callbacks, and the values are the functions that will be
    called at each checkpoint. The functions must receive the wealth distribution
    as an argument.
# Returns
    w_t::Union{Nothing, Matrix{Real}}: Time series of the wealth distribution. Each row is a checkpoint.
    callback_results::Union{Nothing, Dict{Symbol, AbstractVector}}: Results of the callbacks.
"""
function EYSM_base_callbacks(
    N::Integer,
    W_N::Real,
    chi::Real,
    zeta::Real,
    f::Real,
    steps::Integer,
    seed::Integer;
    w0::Union{Nothing, Vector{<:Real}}=nothing,
    initial_conditions::String="uniform",
    save_every::Union{Nothing, Integer}=nothing,
    callbacks::Union{Nothing, Dict{Symbol, Function}}=nothing,
    callbacks_only::Bool=false
)
    # Set the seed
    Random.seed!(seed)
    # Initialize wealth distribution
    w = mc_set_initial_conditions(N, W_N, initial_conditions, w0)
    # Calculate the total wealth
    W = sum(w)
    # Calculate some useful constants
    chif_N = chi * f / N
    zeta_W = zeta / W

    # Initialize the time series
    ## Check the save_every argument
    if isnothing(save_every)
        save_every = N
    end
    ## Check for callbacks
    if callbacks_only
        w_t = nothing
    else
        ## Initialize the time series
        ### Each row of w_t is a checkpoint
        w_t = zeros(typeof(W_N), 1 +(steps ÷ save_every), N)
    end
    # Initialize the callbacks
    callback_results = Dict{Symbol, AbstractVector}()
    if !isnothing(callbacks)
        for (name, f) in callbacks
            callback_results[name] = Vector{typeof(f(w))}(undef, (steps ÷ save_every) + 1)
        end
    end
    # Define a function for saving the results of the callbacks
    function save_callbacks!(w, callbacks, callback_results, idx)
        if !isnothing(callbacks)
            for (name, f) in callbacks
                callback_results[name][idx] = f(w)
            end
        end
    end
    # Define a function for saving the wealth distribution
    function save_wealth!(w_t, w, idx, callbacks_only)
        if !callbacks_only
            @. w_t[idx, :] = w
        end
    end
    # Now save the initial conditions
    idx = 1
    save_callbacks!(w, callbacks, callback_results, idx)
    save_wealth!(w_t, w, idx, callbacks_only)
    # After saving, move to the next index
    idx += 1
    # Simulation loop
    for t in 1:steps
        for exch in 1:N
            # Select two agents
            i, j = rand(1:N, 2)
            # Avoid self-exchanges
            while i == j
                j = rand(1:N)
            end
            # Calculate the wealth exchange
            wi, wj = w[i], w[j]
            δw = Δw(f, wi, wj, zeta_W)
            # Calculate the redistribution terms
            redist_i = EYSM_base_redistribution(wi, W_N, chif_N)
            redist_j = EYSM_base_redistribution(wj, W_N, chif_N)
            # Update the wealth
            w[i] += δw + redist_i
            w[j] += - δw + redist_j
        end
        # Save the wealth distribution and apply the callbacks
        if t % save_every == 0
            # Re normalize the wealth
            w ./= sum(w)
            w .*= W
            save_callbacks!(w, callbacks, callback_results, idx)
            save_wealth!(w_t, w, idx, callbacks_only)
            # After saving, move to the next index
            idx += 1
            # Check for negative wealth
            if any(w .< 0)
                throw(ErrorException("Negative wealth detected. Simulation stopped."))
            end
        end
    end
    # Return the results according to the arguments
    if callbacks_only
        return callback_results
    else
        return w_t, callback_results
    end
end

"""
    EYSM_net_callbacks(g::SimpleGraph,W_N, interaction_mode, taxation_mode, chi, zeta, f,
    steps, seed; w=nothing,initial_conditions="uniform",save_every=nothing,
    callbacks=nothing,callbacks_only=false)
Runs a Monte Carlo simulation of the Extender Yard-Sale model on an arbitrary network and
returns the whole time series of the wealth distribution and the results of the callbacks.

# Arguments
    g::SimpleGraph{<:Integer}: Graph representing the network of agents.
    W_N::Real: Mean wealth of the agents. It is the total wealth divided by the number
    of agents.
    interaction_mode::String: Interaction mode. Options are "A" and "B".
    taxation_mode::String: Taxation mode. Options are "A" and "B".
    chi::Real: Taxation rate.
    zeta::Real: Wealth-Attained-Advantage parameter.
    f::Real: Fraction of the wealth that is redistributed.
    steps::Integer: Number of steps of the simulation, measured in Monte Carlo steps.
    seed::Integer: Seed for the random number generator.
# Optional arguments
    w0::Union{Nothing, Vector{<:Real}}=nothing: Initial wealth distribution.
    Only used if initial_conditions="custom". Default is nothing.
    initial_conditions::String=nothing: Initial condition. Options are
    "uniform", "random", "noisy" and "custom". Default is "uniform". If "custom" is
    chosen, the w argument must be provided.
    save_every::Union{Nothing, Integer}=nothing: Save the wealth distribution every
    save_every steps. Default is nothing, which means saving every N steps.
    callbacks::Union{Nothing, Dict{Symbol, Function}}=nothing: Dictionary of callbacks.
    The keys are the names of the callbacks, and the values are the functions that will be
    called at each checkpoint. The functions must receive the wealth distribution as an
    argument.
# Returns
    w_t::Union{Nothing, Matrix{Real}}: Time series of the wealth distribution.
    Each row is a checkpoint.
    callback_results::Union{Nothing, Dict{Symbol, AbstractVector}}: Results of the callbacks.
"""
function EYSM_net_callbacks(
    g::SimpleGraph{<:Integer},
    W_N::Real,
    interaction_mode::String,
    taxation_mode::String,
    chi::Real,
    zeta::Real,
    f::Real,
    steps::Integer,
    seed::Integer;
    w0::Union{Nothing, Vector{<:Real}}=nothing,
    initial_conditions::String="uniform",
    save_every::Union{Nothing, Integer}=nothing,
    callbacks::Union{Nothing, Dict{Symbol, Function}}=nothing,
    callbacks_only::Bool=false
)

    # Get data from the graph
    N, _, edgelist, gc = get_graph_data(g)
    # Set the seed
    Random.seed!(seed)
    # Initialize wealth distribution
    w = mc_set_initial_conditions(N, W_N, initial_conditions, w0)
    # Calculate the total wealth
    W = sum(w)
    # Calculate some useful constants
    chif_N = chi * f / N
    chif_N2 = chi * f / (N^2)
    zeta_W = zeta / W

    # Initialize the time series
    ## Check the save_every argument
    if isnothing(save_every)
        save_every = N
    end
    ## Check for callbacks
    if callbacks_only
        w_t = nothing
    else
        ## Initialize the time series
        ### Each row of w_t is a checkpoint
        w_t = zeros(typeof(W_N), 1 +(steps ÷ save_every), N)
    end
    # Initialize the callbacks
    callback_results = Dict{Symbol, AbstractVector}()
    if !isnothing(callbacks)
        for (name, f) in callbacks
            callback_results[name] = Vector{typeof(f(w))}(undef, (steps ÷ save_every) + 1)
        end
    end
    # Define a function for saving the results of the callbacks
    function save_callbacks!(w, callbacks, callback_results, idx)
        if !isnothing(callbacks)
            for (name, f) in callbacks
                callback_results[name][idx] = f(w)
            end
        end
    end
    # Define a function for saving the wealth distribution
    function save_wealth!(w_t, w, idx, callbacks_only)
        if !callbacks_only
            @. w_t[idx, :] = w
        end
    end
    # Now save the initial conditions
    idx = 1
    save_callbacks!(w, callbacks, callback_results, idx)
    save_wealth!(w_t, w, idx, callbacks_only)
    # After saving, move to the next index
    idx += 1

    # Best ways to access the graph:
    # Random nodes: Instantiate 1:N
    nodes = 1:N
    # Neighbors of node: Neighbors list or neighbors(g, node) from Graphs.jl
    nl = [neighbors(gc, i) for i in nodes]

    # Check
    if interaction_mode ∉ ["A","B"]
        throw(ArgumentError("Invalid interaction mode. Choose 'A' or 'B'."))
    elseif taxation_mode ∉ ["A","B"]
        throw(ArgumentError("Invalid taxation mode. Choose 'A' or 'B'."))
    end

    # Simulation loop
    # IMA: Randomly select a link
    if interaction_mode == "A"
        # TMA: Tax exchanging agents
        if taxation_mode == "A"
            for t in 1:steps
                for exch in nodes
                    i,j = rand(edgelist)
                    wi, wj = w[i], w[j]
                    # Calculate the wealth exchange
                    δw = Δw(f, wi, wj, zeta_W)
                    # Update the wealth and apply the taxation
                    w[i] += δw - chif_N * wi
                    w[j] += -δw - chif_N * wj
                    # Everybody receives the benefits of the taxation
                    # Broadcasting in TMA is more efficient! (I measured it)
                    @. w += chif_N2 * (wi + wj)
                end
                # Save the wealth distribution
                if t % save_every == 0
                    # Re normalize the wealth
                    w ./= sum(w)
                    w .*= W

                    save_callbacks!(w, callbacks, callback_results, idx)
                    save_wealth!(w_t, w, idx, callbacks_only)
                    # After saving, move to the next index
                    idx += 1
                    # Check for negative wealth
                    if any(w .< 0)
                        throw(ArgumentError("Negative wealth detected. Simulation stopped."))
                        return w_t
                    end
                end
            end
        # TMB: Tax two random agents
        elseif taxation_mode == "B"
            for t in 1:steps
                taxed_agents = [(rand(nodes), rand(nodes)) for _ in 1:N]
                for exch in nodes
                    # Exchange agents
                    i,j = rand(edgelist)
                    wi, wj = w[i], w[j]
                    # Taxed agents
                    i_taxed, j_taxed = taxed_agents[exch]
                    wi_taxed, wj_taxed = w[i_taxed], w[j_taxed]
                    # Calculate the wealth exchange
                    δw = Δw(f, wi, wj, zeta_W)
                    # Update the wealths
                    w[i] += δw
                    w[j] -= δw
                    # Tax the agents
                    w[i_taxed] -= chif_N * wi_taxed
                    w[j_taxed] -= chif_N * wj_taxed
                    # Everybody receives the benefits of the taxation
                    # Trying to reduce the number of memory allocations by using a for loop
                    # It has little impact on performance on TMB (I measured it)
                    tax_benefit = chif_N2 * (wi_taxed + wj_taxed)
                    for i in nodes
                        w[i] += tax_benefit
                    end
                end
                # Save the wealth distribution
                if t % save_every == 0
                    # Re normalize the wealth
                    w ./= sum(w)
                    w .*= W
                    save_callbacks!(w, callbacks, callback_results, idx)
                    save_wealth!(w_t, w, idx, callbacks_only)
                    # After saving, move to the next index
                    idx += 1
                    # Check for negative wealth
                    if any(w .< 0)
                        throw(ArgumentError("Negative wealth detected. Simulation stopped."))
                        return w_t
                    end
                end
            end
        end
    # IMB: Randomly select a node and one of its neighbors
    elseif interaction_mode == "B"
        # TMA: Tax exchanging agents
        if taxation_mode == "A"
            for t in 1:steps
                rand_nodes = rand(1:N, N)
                rand_neighbors = [rand(nl[i]) for i in rand_nodes]
                for exch in nodes
                    i = rand_nodes[exch]
                    j = rand_neighbors[exch]
                    wi, wj = w[i], w[j]
                    # Calculate the wealth exchange
                    δw = Δw(f, wi, wj, zeta_W)
                    # Update the wealth
                    w[i] += δw - chif_N * wi
                    w[j] += -δw - chif_N * wj
                    # Everybody receives the benefits of the taxation
                    # Broadcasting in TMA is more efficient! (I measured it)
                    @. w += chif_N2 * (wi + wj)
                end
                # Save the wealth distribution
                if t % save_every == 0
                    # Re normalize the wealth
                    w ./= sum(w)
                    w .*= W
                    save_callbacks!(w, callbacks, callback_results, idx)
                    save_wealth!(w_t, w, idx, callbacks_only)
                    # After saving, move to the next index
                    idx += 1
                    # Check for negative wealth
                    if any(w .< 0)
                        throw(ArgumentError("Negative wealth detected. Simulation stopped."))
                        return w_t
                    end
                end
            end
        # TMB: Tax two random agents
        elseif taxation_mode == "B"
            for t in 1:steps
                taxed_agents = [(rand(nodes), rand(nodes)) for _ in 1:N]
                rand_nodes = rand(1:N, N)
                rand_neighbors = [rand(nl[i]) for i in rand_nodes]
                for exch in nodes
                    # Exchange agents
                    i = rand_nodes[exch]
                    j = rand_neighbors[exch]
                    wi, wj = w[i], w[j]
                    # Taxed agents
                    i_taxed, j_taxed = taxed_agents[exch]
                    wi_taxed, wj_taxed = w[i_taxed], w[j_taxed]
                    # Calculate the wealth exchange
                    δw = Δw(f, wi, wj, zeta_W)
                    # Update the wealth
                    w[i] += δw
                    w[j] -= δw
                    # Tax the agents
                    w[i_taxed] -= chif_N * wi_taxed
                    w[j_taxed] -= chif_N * wj_taxed
                    # Everybody receives the benefits of the taxation
                    # Trying to reduce the number of memory allocations by using a for loop
                    # It has a little impact on performance on TMB (I measured it)
                    tax_benefit = chif_N2 * (wi_taxed + wj_taxed)
                    for i in nodes
                        w[i] += tax_benefit
                    end
                end
                # Save the wealth distribution
                if t % save_every == 0
                    # Re normalize the wealth
                    w ./= sum(w)
                    w .*= W

                    save_callbacks!(w, callbacks, callback_results, idx)
                    save_wealth!(w_t, w, idx, callbacks_only)
                    # After saving, move to the next index
                    idx += 1
                    # Check for negative wealth
                    if any(w .< 0)
                        throw(ArgumentError("Negative wealth detected. Simulation stopped."))
                        return w_t
                    end
                end
            end
        end
    end
    # Return the results according to the arguments
    if callbacks_only
        return callback_results
    else
        return w_t, callback_results
    end
end
