# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Connectivity{PL,N}

A connectivity list of `N` indices representing a polytope
of type `PL`. Indices are taken from a global vector of
[`Point`](@ref) items.

Connectivity objects are constructed with the [`connect`](@ref)
helper function.
"""
struct Connectivity{PL<:Polytope,N}
  list::NTuple{N,Int}
end

"""
    connect(list, PL)

Connect a `list` of indices from a global vector of points
into a polytope of type `PL`.

## Example
```
Δ = connect((1,2,3), Triangle)
```
"""
connect(list::NTuple{N,Int}, PL::Type{<:Polytope}) where {N} = Connectivity{PL,N}(list)

function materialize(connec::Connectivity{PL},
                     points::AbstractVector{P}) where {PL<:Polytope,P<:Point}
  PL(view(points, collect(connec.list)))
end
