# This file contains some useful functions used to analyze the ODE data
"""
    rescale_t(t::AbstractVector,N::Integer)
Rescale the ODE time vector to MC steps.
# Arguments
    t::AbstractVector: Time vector of the ODE solution.
    N::Integer: Number of nodes in the network.
# Returns
    rescaled_t::AbstractVector: Rescaled time vector.
"""
rescale_t(t::AbstractVector,N::Integer) = 2*(N^2) * t


"""
    get_lambda(g::SimpleGraph{<:Integer}, interaction_mode::String,
    taxation_mode::String, T::Real, x_ss::Vector{<:Real})
Get the Lambda matrix for the ODE system.
# Arguments
    g::SimpleGraph{<:Integer}: Graph object.
    interaction_mode::String: Interaction mode.
    taxation_mode::String: Taxation mode.
    T::Real: Taxation rate.
    x_ss::Vector{<:Real}: Steady state vector.
# Returns
    Lambda::Matrix{Float64}: Lambda matrix.
# Details
    Coming soon.
"""
function get_lambda(
    g::SimpleGraph{<:Integer},
    interaction_mode::String,
    taxation_mode::String,
    T::Real,
    x_ss::Vector{<:Real}
    )
    N = nv(g)
    kappa,beta = get_kappa_beta(g, interaction_mode, taxation_mode)

    # From kappa, vector, to matrix
    kappa_mat = zeros(N,N)
    edgelist = [(e.src,e.dst) for e in edges(g)]
    for (l, (i,j)) in enumerate(edgelist)
        kappa_mat[i,j] = kappa[l]
    end
    # Build the I matrix
    I = zeros(N,N)
    for i in 1:N
        xi = x_ss[i]
        for j in i:N
            xj = x_ss[j]
            b = Float64(xi < xj)
            xj = b
            I[i,j] =  1 - b
        end
    end

    # Build the Lambda matrix
    Lambda = zeros(N,N)
    for i in 1:N
        xi = x_ss[i]
        for j in 1:N
            if i == j
                exch = 0.0
                for k in 1:N
                    xk = x_ss[k]
                    i_ik = I[i,k]
                    exch += kappa_mat[i,k] * ( (i_ik * (2*xi - xk)) +  ( (1 - i_ik) * xk) )
                end

                beta_i = beta[i]
                Lambda[i,j] = beta_i * ( (T/(N^2)) - (T/N) )  + exch
            else
                i_ij = I[i,j]
                Lambda[i,j] = (T/(N^2)) * beta[j] + kappa_mat[i,j] * ( (-i_ij * xi) + ( (1-i_ij) * (xi - 2*x_ss[j]) ) )
            end
        end
    end
    return Lambda
end

"""
    get_max_eigenvalue(Lambda::Array{<:Real})
Get the maximum real part of the eigenvalues of the a real matrix.
We use this function to determine the stability of the ODE system.
# Arguments
    Lambda::Array{<:Real}: Matrix.
# Returns
    max_eigenvalue::Real: Maximum real part of the eigenvalues.
"""
function get_max_eigenvalue(Lambda::Array{<:Real,2})
    return maximum(real(eigen(Lambda).values))
end
