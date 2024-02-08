@computed struct IndexedAdjacenciesTopology{K,dim} <: Topology
    # vertices::AbstractVector{Point{dim, T}}
    simplices::AbstractVector{Connectivity{Ksimplex{K, dim}, K+1}}
    neighboring_simplicies::AbstractVector{NTuple{K+1, Int}}
    R_star_relations::AbstractVector{Int}  # see L De Floriani, A Hui, 2007; Section 5.3

    # the constructor looks a bit crazy because of the @computed macro, see e.g. https://github.com/JuliaLang/julia/issues/18466
    IndexedAdjacenciesTopology(
        simplicies::AbstractVector{C},
        neighboring_simplicies::AbstractVector{NTuple{K_, Int}},
        R_star_relations::AbstractVector{Int}
        ) where{K, K_, dim, C<:Connectivity{Ksimplex{K, dim}}} = new{K,dim}(simplicies, neighboring_simplicies, R_star_relations)
end
