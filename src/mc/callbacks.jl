# Callbacks functions for the Monte Carlo simulations
"""
    get_x1(w::Vector{<:Real})
Calculate the wealth fraction of the richest agent in the system.
# Arguments
    w::Vector{<:Real}: Wealth distribution.
# Returns
    x1::Real: Wealth fraction of the richest agent.
"""
get_x1(w::Vector{<:Real}) = maximum(w) / sum(w)

"""
    get_xi(w::Vector{<:Real}, i::Int)
Calculate the wealth fraction of the i-th richest agent in the system.
# Arguments
    w::Vector{<:Real}: Wealth distribution.
    i::Int: Rank of the agent.
# Returns
    xi::Real: Wealth fraction of the i-th richest agent.
"""
get_xi(w::Vector{<:Real}, i::Int) = sort(w, rev=true)[i] / sum(w)

"""
    get_gini(w::Vector{<:Real})
Calculate the Gini coefficient of the wealth distribution.
# Arguments
    w::Vector{<:Real}: Wealth distribution.
# Details
The Gini coefficient is a measure of statistical dispersion intended to represent the income
or wealth distribution of a nation's residents, and is the most commonly used measure
of inequality. It ranges from 0 (perfect equality) to 1 (perfect inequality).
The formula for the Gini coefficient is:
```math
g = \\frac{1}{2NW} \\sum_{i,j}^{N} |w_i - w_j|
```
where `N` is the number of agents, `W` is the total wealth, and `w_i` is the wealth of
agent `i`.
# Returns
    gini::Real: Gini coefficient.
"""
function get_gini(w::Vector{<:Real})
    g = eltype(w)(0)
    N = length(w)
    for i in 1:N
        wi = w[i]
        for j in 1:N
            g += abs(wi - w[j])
        end
    end
    return g / (2*N*sum(w))
end

"""
    get_lorenz(w::Vector{<:Real})
Calculate the Lorenz curve of the wealth distribution.
# Arguments
    w::Vector{<:Real}: Wealth distribution.
# Details
The Lorenz curve is a graphical representation of the distribution of income or wealth.
It plots the cumulative percentage of total income received against the cumulative
percentage of recipients, starting with the poorest individual or household.
It is used to represent income or wealth inequality. The Lorenz curve is a concave
curve that lies below the diagonal line of perfect equality, which represents a situation
where everyone has the same income or wealth.
# Returns
    lorenz::Vector{<:Real}: Lorenz curve.
"""
get_lorenz(w::Vector{<:Real}) = cumsum(sort(w)) / sum(w)


"""
    get_R(w::Vector{<:Real}, poverty_line::Real)
Calculate the number of agents with wealth above the poverty line.
# Arguments
    w::Vector{<:Real}: Wealth distribution.
    poverty_line::Real: Poverty line.
# Returns
    R::Int: Number of agents with wealth above the poverty line (mean wealth).
"""
get_R(w::Vector{<:Real}) = sum(w .> mean(w))


"""
    get_r(w::Vector{<:Real}, poverty_line::Real)
Calculate the fraction of agents with wealth above the poverty line.
# Arguments
    w::Vector{<:Real}: Wealth distribution.
    poverty_line::Real: Poverty line.
# Returns
    r::Real: Fraction of agents with wealth above the poverty line.
"""
get_r(w::Vector{<:Real}) = typeof(w[1])(get_R(w) / length(w))

"""
    get_u(w::Vector{<:Real}, poverty_line::Real)
Calculate the relative mean wealth of agents with wealth above the poverty line (mean wealth)
# Arguments
    w::Vector{<:Real}: Wealth distribution.
    poverty_line::Real: Poverty line.
# Returns
    u::Real: Mean wealth of agents with wealth above the poverty line.
"""
function get_u(w::Vector{<:Real})
    if isempty(w[w .> mean(w)])
        return typeof(w[1])(0)
    else
        return mean(w[w .> mean(w)])/ sum(w)
    end
end
