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
include("utils/ode_auxiliary_functions.jl")
include("utils/graph_functions.jl")

# ODE
include("ode/ode_solvers.jl")

# MC
include("mc/mc_simulator.jl")
include("mc/callbacks.jl")

# Data analysis
include("data_analysis/ode_analysis_functions.jl")
include("data_analysis/mc_analysis_functions.jl")

# Exported functions
# Utils
export get_kappa_beta
export mc_set_initial_conditions, Δw, ηij, EYSM_base_redistribution
export ode_set_initial_conditions
export get_giant_component, get_graph_data
# ODE
export dxdt_net!, solve_ode_net, solve_ode_net_SS
# MC
export EYSM_base_full, EYSM_net_full
export get_x1, get_xi, get_gini, get_lorenz, get_R, get_u
## Development functions
export EYSM_base_callbacks, EYSM_net_callbacks
# Data analysis
export rescale_t
export get_x1, get_avg_x1

end
