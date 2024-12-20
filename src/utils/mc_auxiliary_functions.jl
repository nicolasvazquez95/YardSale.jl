# mc_auxiliary_functions.jl
# A set of functions to be used in the Monte Carlo simulations.
"""
    mc_set_initial_conditions(N::Integer, W_N::Real, initial_conditions::String,
        w::Union{Nothing, Vector{<:Real}}=nothing)
Set the initial conditions for the Monte Carlo simulation.
# Arguments
    N::Integer: Number of agents.
    W_N::Real: Total wealth per agent.
    initial_conditions::String: Initial condition. Options are "uniform", "random", "noisy"
    and "custom".
    w::Union{Nothing, Vector{<:Real}}=nothing: Initial wealth distribution. Only used if
    initial_conditions="custom". Default is nothing.
# Returns
    w::Vector{Real}: Initial wealth distribution.
"""
function mc_set_initial_conditions(
    N::Integer,
    W_N::Real,
    initial_conditions::String,
    w::Union{Nothing, Vector{<:Real}}=nothing
)
    # Check if the initial conditions are valid
    valid_initial_conditions = ["uniform", "random", "noisy", "custom"]
    if initial_conditions ∉ valid_initial_conditions
        throw(ArgumentError("Invalid initial condition. Must be one of
        $valid_initial_conditions."))
    end
    W = W_N * N
    if initial_conditions == "uniform"
        w = fill(W_N, N)
    elseif initial_conditions == "random"
        w = rand(typeof(W_N),N)
    elseif initial_conditions == "noisy"
        gauss = Normal(0.0, 0.01)
        noise = rand(gauss, N)
        w = fill(W_N, N) + typeof.(W_N).(noise)
    elseif initial_conditions == "custom"
        if isnothing(w)
            throw(ArgumentError("If initial_conditions is set to 'custom',
             w must be provided."))
        end
        if length(w) != N
            throw(ArgumentError("The length of w must be equal to N."))
        end
        # Every element of w must be non-negative
        if any(w .< 0)
            throw(ArgumentError("Every element of w must be non-negative."))
        end
        # Throw a warning if the sum of w is too different from W
        if (sum(w) ≈ W) == false
            @warn "The sum of w is different from W. It will be normalized to W by default."
        end
    end
    # Normalize w
    w .*= W/sum(w)
    return w
end

"""
    Δw(i, j, wi, wj, f)
Calculate the amount of wealth exchanged between agents i and j in the EYSM model.
# Arguments
    i::Integer: Index of agent i.
    j::Integer: Index of agent j.
    wi::Real: Wealth of agent i.
    wj::Real: Wealth of agent j.
    f::Real: Fraction of the minimum wealth exchanged between agents i and j.
# Returns
    Δw::Real: Amount of wealth exchanged between agents i and j.
"""
Δw(i::Integer, j::Integer, wi::Real, wj::Real) = f * min(wi,wj) * ηij(wi, wj, zeta_W)

"""
    ηij(wi, wj, zeta_W)
Calculate the stochastic variable ηij in the EYSM model.
# Arguments
    wi::Real: Wealth of agent i.
    wj::Real: Wealth of agent j.
    zeta_W::Real: Normalized bias towards the richest agent.
# Details
The expected value of ηij is given by
```math
\\langle\\eta_{ij}\\rangle = \\zeta \\frac{w_i - w_j}{W}
```
This variable is used to determine the direction of the wealth exchange between agents
i and j, and can take values ±1.
# Returns
    ηij::Real Stochastic variable ηij.
"""
function ηij(wi::Real, wj::Real, zeta_W::Real)
    # Random stochastic variable
    r = rand()
    # Calculate the bias
    p = zeta_W * (wi - wj)
    bias = (1 + p) / 2
    # If r<bias, return 1, else return -1
    return 2 * (r < bias) - 1
end

"""
    EYSM_base_redistribution(wi::Real, W_N::Real, chif_N::Real)
Calculates the redistribution term in the EYSM model for an agent i with wealth wi.
# Arguments
    wi::Real: Wealth of agent i.
    W_N::Real: Total wealth per agent.
    chif_N::Real: Taxation and redistribution rate.
# Details
In the EYSM model, the redistribution term of the wealth update of agent i ``is`` given by
```math
\\frac{\\chi f}{N}(\\frac{W}{N} - w_i)
```
where f is the fraction of the minimum wealth exchanged between agents, and ``f``
is playing the role of the Δt term in the Fokker-Planck approach, as the
temporal unit of time of the Monte Carlo simulation. ``chi`` represents the taxation and
redistribution rate and N is the number of agents.
# Returns
    redist::Real: Redistribution term for agent i.
"""
EYSM_base_redistribution(wi::Real, W_N::Real, chif_N::Real) = chif_N * (W_N - wi)
