# kappa_beta.jl
# Here we define the functions to get the κ matrix and the β vector for solving the ODEs.
# Each Interaction Mode and Taxation Mode has its own function to get the κ matrix and the β
# vector.
# The main function get_C_b receives the graph, the Interaction Mode and the Taxation Mode
# and returns the κ matrix and the β vector.

### Private functions
# IMA/TMA
function _get_kappa_beta_ima_tma(N,l,k)
    kappa = zeros(l)
    l1 = 1/l

    @. kappa = l1
    beta = k/l
    return kappa,beta
end

# IMA/TMB
function _get_kappa_beta_ima_tmb(N,l)
    kappa = zeros(l)
    l1 = 1/l

    @. kappa = l1
    beta = 2*ones(N) / N
    return kappa,beta
end

# IMB/TMA
function _get_kappa_beta_imb_tma(N,edgelist,k,A)
    l = size(edgelist, 1)
    kappa = zeros(l)
    beta = zeros(N)

    # Use broadcasting for faster computation
    @. kappa = (1 / k[edgelist[:, 1]]) + (1 / k[edgelist[:, 2]])
    kappa /= N

    # z_i = sum_j A_ij / k_j
    z = vec(sum(A ./ k, dims=1)[:]')
    # beta_i = (1/N) * (1 + z_i)
    @. beta = (1/N) * (1 + z)

    return kappa,beta
end

# IMB/TMB
function _get_kappa_beta_imb_tmb(N,edgelist,k)
    l = size(edgelist, 1)
    kappa = zeros(l)

    # Use broadcasting for faster computation
    @. kappa = (1 / k[edgelist[:, 1]]) + (1 / k[edgelist[:, 2]])
    kappa /= N

    beta = 2*ones(N) / N
    return kappa,beta
end

### Main function (public)
"""
    get_kappa_beta(g::SimpleGraph{<:Integer},interaction_mode::String,taxation_mode::String)
Return the κ matrix and the β vector for the ODEs, given a graph,
an Interaction Mode and a Taxation Mode.
# Arguments
    g::SimpleGraph{<:Integer}: Undirected graph.
    interaction_mode::String: Interaction Mode. Can be "A" or "B".
    taxation_mode::String: Taxation Mode. Can be "A" or "B".
# Returns
    Tuple{Vector{Float64},Vector{Float64}}: κ matrix and β vector.
# Examples
```julia
g = complete_graph(4)
YardSale.get_kappa_beta(g,"A","A")
```
"""
function get_kappa_beta(
    g::SimpleGraph{<:Integer},
    interaction_mode::String,
    taxation_mode::String
    )
    if interaction_mode=="A"
        l = length(edges(g))
        N = nv(g)
        if taxation_mode=="A"
            return _get_kappa_beta_ima_tma(N,l,degree(g))
        elseif taxation_mode=="B"
            return _get_kappa_beta_ima_tmb(N,l)
        else
            throw(ArgumentError("Invalid taxation mode."))
        end
    elseif interaction_mode=="B"
        N = nv(g)
        # Edgelist is a Lx2 matrix where each row is an edge
        edgelist = transpose(hcat([[e.src, e.dst] for e in edges(g)]...))
        k = degree(g)
        if taxation_mode=="A"
            # A is the adjacency matrix A_ij
            A = adjacency_matrix(g)
            return _get_kappa_beta_imb_tma(N,edgelist,k,A)
        elseif taxation_mode=="B"
            return _get_kappa_beta_imb_tmb(N,edgelist,k)
        else
            throw(ArgumentError("Invalid taxation mode."))
        end
    else
        throw(ArgumentError("Invalid interaction mode."))
    end
    return nothing
end
