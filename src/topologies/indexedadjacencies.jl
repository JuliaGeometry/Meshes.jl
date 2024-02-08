# Implement the IndexedAdjacenciesTopologies (IA) described in
# - Paoluzzi, A., Bernardini, F., Cattani, C., & Ferrucci, V. (1993). Dimension-independent modeling with simplicial complexes. ACM Transactions on Graphics (TOG), 12(1), 56-102.
# and compared in
# De Floriani, L., & Hui, A. (2007). Shape Representations Based on Simplicial and Cell Complexes. In Eurographics (State of the Art Reports) (pp. 63-87).
# Assumptions
# - topology is made up of one type of simplex only
# - everything is strongly connected

function IndexedAdjacenciesTopology(connections::AbstractVector{<:Connectivity{Ksimplex{K, dim}, N}}) where {K, dim, N}
    simplicies = connections

    # find neighboring simplicies with max K
    neighboring_simplicies = NTuple{K+1, Int}[]; sizehint!(neighboring_simplicies, length(connections))
    # preallocate
    neighboring_simplicies_for_current_simplex = Int[]; sizehint!(neighboring_simplicies_for_current_simplex, K+1)
    # very inefficient way to find neighbors.
    for (idx_self, simpl_self) in enumerate(simplicies)
        # repeatedly find simplex that shares all but one vertex
        for leftout_vertex in indices(simpl_self)
            joint_vertices = [i for i in indices(simpl_self) if i != leftout_vertex]
            idx = findall(simpl_other->all(joint_vertices .∈ [indices(simpl_other)]),
                          simplicies)
            @assert length(idx) ∈ [1,2] "$(idx)"
            if length(idx) == 2  # self and other
                push!(neighboring_simplicies_for_current_simplex, first(idx[idx .!= idx_self]))
            elseif length(idx) == 1  # no other
                push!(neighboring_simplicies_for_current_simplex, -1)
            end
        end
        push!(neighboring_simplicies, NTuple{K+1}(neighboring_simplicies_for_current_simplex))
        empty!(neighboring_simplicies_for_current_simplex)
    end

    # construct R_star
    all_vertex_indices = unique(reduce(vcat, collect.(indices.(simplicies))))
    R_star_relations=Int[
        findfirst(simpl -> idx ∈ indices(simpl), simplicies)
        for idx in all_vertex_indices
    ]

    IndexedAdjacenciesTopology(simplicies, neighboring_simplicies, R_star_relations)
end

function setuptest()
    pts = [Point(0,0,0), Point(-1,-1,0), Point(-1,1,0), Point(1,1,0), Point(1,-1,0)]
    cons = connect.([(1, 2, 3), (1, 3, 4), (1, 4, 5), (1, 5, 2)], Ksimplex{2, 3})
    topo = IndexedAdjacenciesTopology(cons)
end

paramdim(::IndexedAdjacenciesTopology{K}) where {K} = K

# ---------------------
# HIGH-LEVEL INTERFACE
# ---------------------

nvertices(t::IndexedAdjacenciesTopology) = length(t.R_star_relations)

# function faces(t::SimpleTopology, rank)
#   cs = t.connec
#   (cs[i] for i in 1:length(cs) if paramdim(cs[i]) == rank)
# end

element(t::IndexedAdjacenciesTopology, idx) = t.simplicies[idx]

nelements(t::IndexedAdjacenciesTopology) = length(t.simplicies)
