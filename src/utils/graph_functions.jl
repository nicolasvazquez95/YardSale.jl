# In this script, we define functions to analyze graphs. Mainly, we want functions to get
# the giant component of a graph.
"""
    get_giant_component(g::SimpleGraph)
Returns the giant component of a graph.
"""
function get_giant_component(g::SimpleGraph)
    cc = connected_components(g)
    max_cc = argmax(length, cc)
    return cc[max_cc]
end

"""
    get_graph_data(g::SimpleGraph)
Returns the number of nodes, the number of links,
and the edgelist of the giant component of a graph.

# Details
If the graph has more than one connected component, a warning is issued,
and the giant component is used for the simulation.
"""
function get_graph_data(g::SimpleGraph)
    # Get the giant component
    gc = get_giant_component(g)
    if nv(gc) ≠ nv(g)
        @warn "The graph has a giant component of $(nv(gc)) nodes.
        The original graph has $(nv(g)) nodes.
        We will use the giant component for the simulation."
    end
    # Get the number of nodes
    N = nv(gc)
    # Get the number of links
    l = length(edges(gc))
    # Get the edgelist
    edgelist = [[e.src,e.dst] for e in edges(gc)]
    return N, l, edgelist, gc
end
