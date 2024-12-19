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

# ODE
include("ode/ode_solvers.jl")

# Data analysis
include("data_analysis/analysis_functions.jl")

# Exported functions
# Utils
export get_kappa_beta
# ODE
export dxdt_net!, solve_ode_net, solve_ode_net_steady_state

# Data analysis
export rescale_t

end
