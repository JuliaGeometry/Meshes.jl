# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Connectivity{F,N}

A connectivity list of `N` indices representing a [`Polytope`](@ref)
of type `F`. Indices are taken from a global vector of [`Point`](@ref).

Connectivity objects are constructed with the [`connect`](@ref) function.
"""
struct Connectivity{F<:Polytope,N}
  list::NTuple{N,Int}
end

"""
    facetype(connectivity)

Return the face type (e.g. Triangle) of the `connectivity`.
"""
facetype(::Type{Connectivity{F,N}}) where {F,N} = F
facetype(c::Connectivity) = facetype(typeof(c))

"""
    connect(list, F)

Connect a `list` of indices from a global vector of [`Point`](@ref)
into a [`Polytope`](@ref) of type `F`.

## Example
```
Δ = connect((1,2,3), Triangle)
```
"""
connect(list::NTuple{N,Int}, F::Type{<:Polytope}) where {N} = Connectivity{F,N}(list)

"""
    materialize(connec, points)

Materialize a face using the `connec` list and a global vector of `points`.
"""
function materialize(connec::Connectivity{F},
                     points::AbstractVector{P}) where {F<:Polytope,P<:Point}
  F(view(points, collect(connec.list)))
end

function Base.show(io::IO, c::Connectivity{F}) where {F}
  print(io, "$F$(c.list)")
end
