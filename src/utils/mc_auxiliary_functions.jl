# mc_auxiliary_functions.jl
# A set of functions to be used in the Monte Carlo simulations.
"""
    mc_set_initial_conditions(N, W_N, initial_conditions,w=nothing)
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
    mc_set_beta(N, initial_conditions, beta0=nothing)
Set the initial risk aversion for the Monte Carlo simulation.

# Arguments
    N::Integer: Number of agents.
    initial_conditions::String: Initial condition. Options are "uniform", "random", "noisy"
    and "custom".
    beta0::Union{Nothing, Vector{<:Real}}=nothing: Initial risk propension. Only used if
    initial_conditions="custom". Default is nothing.

# Returns
    beta_::Vector{Real}: Initial risk aversion distribution.
"""
function mc_set_beta(
    N::Integer,
    beta::String,
    beta0::Union{Nothing, Real, Vector{<:Real}}=nothing;
    σ::Real=0.01
)
    # Check if the initial conditions are valid
    valid_initial_conditions = ["uniform", "random", "noisy", "custom"]
    if beta ∉ valid_initial_conditions
        throw(ArgumentError("Invalid initial condition. Must be one of
        $valid_initial_conditions."))
    end
    if beta == "uniform"
        beta_ = fill(beta0, N)
    elseif beta == "random"
        beta_ = rand(N)
    elseif beta == "noisy"
        gauss = Normal(beta0, σ)
        noise = rand(gauss, N)
        beta_ = fill(beta0, N) + noise
    elseif beta == "custom"
        if isnothing(beta0)
            throw(ArgumentError("If initial_conditions is set to 'custom',
             beta0 must be provided."))
        end
        if length(beta0) != N
            throw(ArgumentError("The length of beta0 must be equal to N."))
        end
        # Every element of beta0 must be non-negative
        if any(beta0 .< 0)
            throw(ArgumentError("Every element of beta0 must be non-negative."))
        end
        beta_ = beta0
    end
    return beta_
end


"""
    Δw(f, wi, wj)
Calculate the amount of wealth exchanged between agents i and j in the EYSM model.
# Arguments
    f::Real: Fraction of the minimum wealth exchanged between agents i and j.
    wi::Real: Wealth of agent i.
    wj::Real: Wealth of agent j.

# Returns
    Δw::Real: Amount of wealth exchanged between agents i and j.
"""
Δw(f::Real, wi::Real, wj::Real,zeta_W::Real) = f * min(wi,wj) * ηij(wi, wj, zeta_W)


"""
    Δw(wi,wj,ri,rj)
Calculate the amount of wealth exchanged between agents i and j in the YS model
with risk aversion.
# Arguments
    wi::Real: Wealth of agent i.
    wj::Real: Wealth of agent j.
    ri::Real: Risk propension of agent i.
    rj::Real: Risk propension of agent j.
# Returns
    Δw::Real: Amount of wealth exchanged between agents i and j.
"""
Δw_risk(wi::Real, wj::Real, ri::Real, rj::Real) = min(ri*wi, rj*wj)

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
    EYSM_base_redistribution(wi, W_N, chif_N)
Calculates the redistribution term in the EYSM model for an agent i with wealth wi.
# Arguments
    wi::Real: Wealth of agent i.
    W_N::Real: Total wealth per agent.
    chif_N::Real: Taxation and redistribution rate.
# Details
In the EYSM model, the redistribution term of the wealth update of agent ``i`` is given by
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
