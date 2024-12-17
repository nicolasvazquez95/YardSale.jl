# This test is for the kappa_beta functions in the kappa_beta.jl file.
# The kappa_beta functions are used to get the κ matrix and the β vector for solving
# the ODEs.
# In this test we check if the kappa_beta functions return the correct values for the
# κ matrix and the β vector. We know that for a fully connected graph the κ matrix
# should be a vector of 1/l and the β vector should be a vector of 1/N, for all the
# Interaction Modes and Taxation Modes.

using YardSale, Test, Graphs

# Fully connected graph
N = 16
g = complete_graph(N)

# Expected values of kappa and beta
l = length(edges(g))
# Dumb test to check if the graph is fully connected
@test l == N*(N-1)/2

# Expected values of kappa and beta (see paper)
expected_kappa = [1/l for i in 1:l]
expected_beta = [2/N for i in 1:N]

# Interaction Modes
IM = ["A","B"]
# Taxation Modes
TM = ["A","B"]
im_tm = [(im,tm) for im in IM for tm in TM]
@testset "Interaction/Taxation Modes" begin
    for (im,tm) in im_tm
        kappa,beta = get_kappa_beta(g,im,tm)
        @test kappa ≈ expected_kappa
        @test beta ≈ expected_beta
    end
end

# Error cases
@testset "Error cases" begin
    @test_throws ArgumentError get_kappa_beta(g,"A","C")
    @test_throws ArgumentError get_kappa_beta(g,"C","A")
end
