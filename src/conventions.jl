# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
Convention used to reconstruct the faces of a polytope
from a list of vertices.
"""
abstract type OrderingConvention end

"""
    connectivities(polytope, rank, convention)

Return the list of `Connectivity` elements used to construct all
faces of rank `rank` of a polytope, following an ordering convention.
"""
function connectivities(polytope::Type{<:Polytope}, rank::Integer,
                        convention::OrderingConvention)
  connectivities(polytope, Val(rank), convention)
end

# ---------------------------
# CONVENTION IMPLEMENTATIONS
# ---------------------------

include("conventions/gmsh.jl")

function facets(p::Polytope{N}, convention=GMSH) where {N}
  faces(p, N-1, convention)
end

function faces(p::Polytope, rank, convention=GMSH)
  (materialize(ord, p.vertices) for ord in connectivities(typeof(p), rank, convention))
end
