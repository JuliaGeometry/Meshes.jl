using DataStructures: Queue

# assumptions
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

import Meshes: paramdim
paramdim(::IndexedAdjacenciesTopology{K}) where {K} = K


function (𝒜::Adjacency{K,K,T})(simpl_index::Int) where {K,T<:IndexedAdjacenciesTopology{K}}
    𝒜.topology.neighboring_simplicies[simpl_index]
end

function (𝒞::Coboundary{0,K,K,T})(vert_idx::Int) where {K, T<:IndexedAdjacenciesTopology}
    𝒜 = Adjacency{K}(𝒞.topology)
    results = Int[]
    indices_to_explore = Queue{Int}()
    enqueue!(indices_to_explore, 𝒞.topology.R_star_relations[vert_idx])
    while !isempty(indices_to_explore)
        idx = dequeue!(indices_to_explore)
        if idx == -1 || idx ∈ results
            continue
        elseif vert_idx ∈ indices(𝒞.topology.simplices[idx])
            push!(results, idx)
            for idx in 𝒜(idx)
                enqueue!(indices_to_explore, idx)
            end
        end
    end
    results
end
