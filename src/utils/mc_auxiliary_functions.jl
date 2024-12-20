# mc_auxiliary_functions.jl
# A set of functions to be used in the Monte Carlo simulations.
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
        # Throw a warning if the sum of w is too different from W
        if (sum(w) ≈ W) == false
            @warn "The sum of w is different from W. It will be normalized to W by default."
        end
    end
    # Normalize w
    w .*= W/sum(w)
    return w
end
