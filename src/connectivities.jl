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
    paramdim(connectivity)

Return the parametric dimension of the `connectivity`.
"""
paramdim(::Type{Connectivity{PL,N}}) where {PL,N} = paramdim(PL)
paramdim(c::Connectivity) = paramdim(typeof(c))

"""
    indices(connectivity)

Return the list of indices of the `connectivity`.
"""
indices(c::Connectivity) = c.list

"""
    connect(list, PL)

Connect a `list` of indices from a global vector of [`Point`](@ref)
into a [`Polytope`](@ref) of type `PL`.

## Example
```
Δ = connect((1,2,3), Triangle)
```
"""
connect(list::Tuple, PL::Type{<:Polytope}) = Connectivity{PL,length(list)}(list)

"""
    materialize(connec, points)

Materialize a face using the `connec` list and a global vector of `points`.
"""
function materialize(connec::Connectivity{PL},
                     points::AbstractVector{P}) where {PL<:Polytope,P<:Point}
  PL(view(points, SVector(connec.list...)))
end

function Base.show(io::IO, c::Connectivity{PL}) where {PL}
  name = prettyname(PL)
  print(io, "$name$(c.list)")
end
