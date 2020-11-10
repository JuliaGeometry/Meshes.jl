# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Connectivity{PL,N}

A connectivity list of `N` indices representing a [`Polytope`](@ref)
of type `PL`. Indices are taken from a global vector of [`Point`](@ref).

Connectivity objects are constructed with the [`connect`](@ref) function.
"""
struct Connectivity{PL<:Polytope,N}
  list::NTuple{N,Int}
end

"""
    polytopetype(connectivity)

Return the face type (e.g. Triangle) of the `connectivity`.
"""
polytopetype(::Type{Connectivity{PL,N}}) where {PL,N} = PL
polytopetype(c::Connectivity) = polytopetype(typeof(c))

"""
    connect(list, PL)

Connect a `list` of indices from a global vector of [`Point`](@ref)
into a [`Polytope`](@ref) of type `PL`.

## Example
```
Δ = connect((1,2,3), Triangle)
```
"""
connect(list::NTuple{N,Int}, PL::Type{<:Polytope}) where {N} = Connectivity{PL,N}(list)

"""
    materialize(connec, points)

Materialize a face using the `connec` list and a global vector of `points`.
"""
function materialize(connec::Connectivity{PL},
                     points::AbstractVector{P}) where {PL<:Polytope,P<:Point}
  PL(view(points, collect(connec.list)))
end

function Base.show(io::IO, c::Connectivity{PL}) where {PL}
  print(io, "$PL$(c.list)")
end

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

function facets(p::Polytope{N}, ordering=GMSH) where {N}
  faces(p, N-1, ordering)
end

function faces(p::Polytope{N,Dim}, rank, ordering=GMSH) where {N,Dim}
  (materialize(ord, p.vertices) for ord in connectivity(typeof(p), ordering, rank))
end
