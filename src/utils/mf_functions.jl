# Analytical functions for the mean field theory

"""
    x1_meanField(T::Real)
Mean field prediction for the wealth of the richest agent in the system as a function of the
temperature T. This is the thermodinamic limit case, i.e. N -> ∞.
# Arguments
    `T::Real`: temperature of the system
# Returns
    `x1::Real`: mean field prediction for the wealth of the richest agent
"""
function x1_meanField(T::Real)
    if T > 1.0
        return 0.0
    else
        return 1.0 - T
    end
end

"""
    x1_meanField(T::Real,N::Int)
Mean field prediction for the wealth of the richest agent in the system as a function of the
temperature T. This is the finite size case, i.e. N is finite.
# Arguments
    `T::Real`: temperature of the system
    `N::Int`: number of agents in the system
# Returns
    `x1::Real`: mean field prediction for the wealth of the richest agent
"""
function x1_meanField(T::Real,N::Int)
    if T > 1.0
        return 1/N
    else
        return 1.0 - ( ((N-1)/N)^2 * T )
    end
end
