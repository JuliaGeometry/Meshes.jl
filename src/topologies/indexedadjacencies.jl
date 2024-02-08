"""
    IndexedAdjacenciesTopology{K, dim, K_} <: Topology

IndexedAdjacenciesTopologies (IA) described in
- Paoluzzi, A., Bernardini, F., Cattani, C., & Ferrucci, V. (1993). Dimension-independent modeling with simplicial complexes. ACM Transactions on Graphics (TOG), 12(1), 56-102.
and compared in
De Floriani, L., & Hui, A. (2007). Shape Representations Based on Simplicial and Cell Complexes. In Eurographics (State of the Art Reports) (pp. 63-87).
Assumptions
- topology is made up of one type of simplex only
- everything is strongly connected

The data structure essentially stores all simplicies of the highest rank (d) and also stores the
neighbors for each simplex that shares a (d-1) face.
For example, if you have triangles (K=2), they will share a segment (K=1) with another triangle,
and each triangle will have exactly (3==K+1) neighbors (although the neighbors may be empty,
which is denoted by the index -1).
Notice that this is all independent of the dimension, as long as K <= dim.

Further, we have the convenience type K_==(K+1) always.
We tried using the @computed macro from ComputedFieldTypes.jl but it's a bit messy.
"""
struct IndexedAdjacenciesTopology{K, dim, K_} <: Topology
    # vertices::AbstractVector{Point{dim, T}}
    simplicies::AbstractVector{Connectivity{Ksimplex{K, dim}, K_}}
    neighboring_simplicies::AbstractVector{NTuple{K_, Int}}
    R_star_relations::AbstractVector{Int}  # see L De Floriani, A Hui, 2007; Section 5.3

    """
    K_ must be equal to K+1
    We tried using the @computed macro from ComputedFieldTypes.jl but it's a bit messy.
    """
    function IndexedAdjacenciesTopology(
        simplicies::AbstractVector{C},
        neighboring_simplicies::AbstractVector{NTuple{K_, Int}},
        R_star_relations::AbstractVector{Int}
        ) where{K, K_, dim, C<:Connectivity{Ksimplex{K, dim}}}

        @assert K+1 == K_
        new{K,dim,K+1}(simplicies, neighboring_simplicies, R_star_relations)
    end
end

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
