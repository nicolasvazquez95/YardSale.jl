"""
    YardSale.jl
`YardSale.jl` is an extensible package for the study of the dynamics of wealth exchange
models in complex networks.
"""
module YardSale
# Packages
using DifferentialEquations
using Graphs
using Distributions
using Random
using Statistics

# Include statements
# Utils
include("utils/kappa_beta.jl")
include("utils/mc_auxiliary_functions.jl")

# ODE
include("ode/ode_solvers.jl")

# Data analysis
include("data_analysis/analysis_functions.jl")

# Exported functions
# Utils
export get_kappa_beta
export mc_set_initial_conditions
# ODE
export dxdt_net!, solve_ode_net, solve_ode_steady_state

# Data analysis
export rescale_t

end
