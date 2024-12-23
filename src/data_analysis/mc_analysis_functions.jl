# This file contains some useful functions used to analyze the MC data
"""
    get_x1(x::AbstractMatrix, dims::Union{Integer, Tuple{Integer, Integer}}=2)
    Get the richest node at each time step.
# Arguments
    x::AbstractMatrix: Matrix of the MC solution.
    dims::Union{Integer, Tuple{Integer, Integer}}: Dimension(s) to get the richest node.
    It can be an integer or a tuple of integers.
# Returns
    x1::AbstractVector: Richest node at each time step.
"""
get_x1(x::AbstractMatrix, dims::Vararg{Integer}=2) = maximum(x, dims=dims)
