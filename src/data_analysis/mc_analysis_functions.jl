# This file contains some useful functions used to analyze the MC data
"""
    get_x1(x::AbstractMatrix, dims::Union{Integer, Tuple{Integer, Integer}}=2)
    Get the richest node at each time step.
# Arguments
    x::AbstractMatrix: Matrix of the MC solution.
    dims::Union{Integer, Tuple{Integer, Integer}}: Dimension(s) to get the richest node.
    It can be an integer or a tuple of integers.
# Returns
    x1::AbstractArray: Richest node at each time step.
"""
function get_x1(x::AbstractArray, dims::Vararg{Integer}=2)
    return maximum(x, dims=dims)
end


"""
    get_avg_x1(x::AbstractMatrix, dims::Vararg{Integer}=2)
Get the temporal average of the richest node and the variance.

# Arguments
    x1::AbstractMatrix: Matrix of the MC simulation.
    dims::Vararg{Integer}: Dimension(s) to get the richest node.
    It can be an integer or a tuple of integers.
# Returns
    avg_x1::Tuple{Float64, Float64}: Tuple with the average and variance of the richest node.
"""
function get_avg_x1(x1::AbstractArray, dims::Integer)
    return mean(x1, dims=dims), var(x1, dims=dims)
end
