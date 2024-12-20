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
w_i \\leftarrow w_i - \\Delta w + \\chi (\\frac{W}{N} - w_i)
w_j \\leftarrow w_j + \\Delta w + \\chi (\\frac{W}{N} - w_j)
````
where ``\\chi`` represents the taxation and redistribution rate.

In the original paper, the authors work analitically, using a Fokker-Planck approach. In this
implementation, we use a Monte Carlo simulation to obtain similar results.

# Examples
    TODO
"""
function EYSM_FCN_full(
    N::Integer,
    W_N::Real,
    chi::Real,
    zeta::Real,
    f::Real,
    steps::Integer,
    seed::Integer;
    w::Union{Nothing, Vector{<:Real}}=nothing,
    initial_conditions::String="uniform",
    save_every::Union{Nothing, Integer}=nothing
)

    # Set the seed
    Random.seed!(seed)

    # 20/12 TODO: We're stuck with the mc_set_initial_conditions function. Implement and test.
    # w = mc_set_initial_conditions(N, W_N, initial_conditions, w)
    # Then, proceed with the rest of the function.
end
