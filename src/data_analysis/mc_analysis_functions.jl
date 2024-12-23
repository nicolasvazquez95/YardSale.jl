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
get_x1(x::AbstractMatrix, dims::Union{Integer, Tuple{Integer, Integer}}=2) = maximum(x, dims=dims)

"""
    get_avg_x1(x::AbstractMatrix, start::Integer, dims::Union{Integer, Tuple{Integer, Integer}}=2)
    Get the time average and variance of the richest node after a given time step.
# Arguments
    x::AbstractMatrix: Matrix of the MC solution.
    start::Integer: Time step to start the average.
    dims::Union{Integer, Tuple{Integer, Integer}}: Dimension(s) to get the richest node.
    It can be an integer or a tuple of integers.
# Returns
    avg_x1::Float64: Time average of the richest node after a given time step.
    var_x1::Float64: Variance of the richest node after a given time step.
"""
function get_avg_x1(x::AbstractMatrix, start::Integer, dims::Union{Integer, Tuple{Integer, Integer}}=2)
    x1 = get_x1(x, dims)
    return mean(x1[start:end]), var(x1[start:end])
end
