# This file contains some useful functions used to analyze the ODE data
"""
    rescale_t(t::AbstractVector,N::Integer)
Rescale the ODE time vector to MC steps.
# Arguments
    t::AbstractVector: Time vector of the ODE solution.
    N::Integer: Number of nodes in the network.
# Returns
    rescaled_t::AbstractVector: Rescaled time vector.
"""
rescale_t(t::AbstractVector,N::Integer) = 2*(N^2) * t
