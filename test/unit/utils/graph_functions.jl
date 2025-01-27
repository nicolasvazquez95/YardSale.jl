# Test for graph_functions.jl
using YardSale, Test, Graphs

@testset "giant_component" begin
    # Test for the get_giant_component function
    # Create a graph
    # We know this is a disconnected graph
    g = erdos_renyi(128, 8/127, seed=10)
    # Get the giant component
    gc = get_giant_component(g)
    # Test that the giant component is a SimpleGraph
    @test isa(gc, SimpleGraph)
end
