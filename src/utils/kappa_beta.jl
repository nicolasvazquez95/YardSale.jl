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
    get_kappa_beta(g::SimpleGraph{<:Integer},IM::String,TM::String)
Return the κ matrix and the β vector for the ODEs, given a graph,
an Interaction Mode and a Taxation Mode.
# Arguments
    g::SimpleGraph{<:Integer}: Undirected graph.
    IM::String: Interaction Mode. Can be "A" or "B".
    TM::String: Taxation Mode. Can be "A" or "B".
# Example
```julia
g = complete_graph(4)
YardSale.get_kappa_beta(g,"A","A")
```
"""
function get_kappa_beta(g::SimpleGraph{<:Integer},IM::String,TM::String)
    if IM=="A"
        l = length(edges(g))
        N = nv(g)
        if TM=="A"
            return _get_kappa_beta_ima_tma(N,l,degree(g))
        elseif TM=="B"
            return _get_kappa_beta_ima_tmb(N,l)
        else
            throw(ArgumentError("Invalid TM"))
        end
    elseif IM=="B"
        N = nv(g)
        # Edgelist is a Lx2 matrix where each row is an edge
        edgelist = transpose(hcat([[e.src, e.dst] for e in edges(g)]...))
        k = degree(g)
        if TM=="A"
            # A is the adjacency matrix A_ij
            A = adjacency_matrix(g)
            return _get_kappa_beta_imb_tma(N,edgelist,k,A)
        elseif TM=="B"
            return _get_kappa_beta_imb_tmb(N,edgelist,k)
        else
            throw(ArgumentError("Invalid TM"))
        end
    else
        throw(ArgumentError("Invalid IM"))
    end
    return nothing
end
