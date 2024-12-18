"""
    EYSM_FCN_full(N::Integer, W_N::Real, chi::Real, zeta::Real, f::Real, steps::Integer,
    seed::Integer; w::Union{Nothing, Vector{<:Real}}=nothing,
    initial_conditions::String="uniform",
    save_every::Union{Nothing, Integer}=nothing)
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
    w::Union{Nothing, Vector{<:Real}}=nothing: Initial wealth distribution.
    Only used if initial_conditions="custom". Default is nothing.
    initial_conditions::String=nothing: Initial condition. Options are
    "uniform", "random", "noisy" and "custom". Default is "uniform". If "custom" is
    chosen, the w argument must be provided.
    save_every::Union{Nothing, Integer}=nothing: Save the wealth distribution every
    save_every steps. Default is nothing, which means saving every N steps.
# Returns

# Details

# Examples
"""
