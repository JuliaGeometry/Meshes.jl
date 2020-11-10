# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
Convention used to reconstruct the faces of a polytope
from a list of vertices.
"""
abstract type OrderingConvention end

"""
    connectivity(polytope, convention, k)

Return the list of `Connectivity` elements used to construct
all `k`-faces of a polytope, following an ordering convention.
"""
function connectivity(polytope::Type{<:Polytope}, convention::OrderingConvention, k::Integer)
  connectivity(polytope, convention, Val(k))
end

# ---------------------------
# CONVENTION IMPLEMENTATIONS
# ---------------------------

include("conventions/gmsh.jl")

function facets(p::Polytope{N}, ordering=GMSH) where {N}
  faces(p, N-1, ordering)
end

function faces(p::Polytope{N,Dim}, rank, ordering=GMSH) where {N,Dim}
  (materialize(ord, p.vertices) for ord in connectivity(typeof(p), ordering, rank))
end
