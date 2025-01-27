# mc_auxiliary_functions.jl
# A set of functions to be used in the ODE simulations.
"""
    ode_set_initial_conditions(N, initial_conditions,x0=nothing)
Set the initial conditions for the ODE simulation.
# Arguments
    N::Integer: Number of agents.
    initial_conditions::String: Initial condition. Options are "uniform", "random", "noisy"
    and "custom".
    x0::Union{Nothing, Vector{<:Real}}=nothing: Initial wealth distribution. Only used if
    initial_conditions="custom". Default is nothing.
# Returns
    x::Vector{Real}: Initial wealth distribution.
"""
function ode_set_initial_conditions(
    N::Integer,
    initial_conditions::String,
    x0::Union{Nothing, Vector{<:Real}}=nothing
)
    # Check if the initial conditions are valid
    valid_initial_conditions = ["uniform", "random", "noisy", "custom"]
    if initial_conditions ∉ valid_initial_conditions
        throw(ArgumentError("Invalid initial condition. Must be one of
        $valid_initial_conditions."))
    end
    W = 1.0
    if initial_conditions == "uniform"
        x = fill(1/N, N)
    elseif initial_conditions == "random"
        x = rand(N)
    elseif initial_conditions == "noisy"
        gauss = Normal(0.0, 0.01)
        noise = rand(gauss, N)
        x = (1/N) * (ones(N) + noise)
    elseif initial_conditions == "custom"
        if isnothing(x0)
            throw(ArgumentError("If initial_conditions is set to 'custom',
             x0 must be provided."))
        end
        if length(x0) != N
            throw(ArgumentError("The length of x0 must be equal to N."))
        end
        # Every element of w must be non-negative
        if any(x0 .< 0)
            throw(ArgumentError("Every element of x0 must be non-negative."))
        end
        # Throw a warning if the sum of w is too different from W
        if (sum(x0) ≈ W) == false
            @warn "The sum of x0 is different from 1. It will be normalized by default."
        end
        x = x0
    end
    # Normalize x
    x ./= W
    return x
end
