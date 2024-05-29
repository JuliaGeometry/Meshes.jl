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
  indices::NTuple{N,Int}

  function Connectivity{PL,N}(indices) where {PL,N}
    assertion(nvertices(PL) == N, lazy"cannot create a $PL with $N vertices")
    if PL <: Ngon
      assertion(N ≥ 3, "Ngon requires 3 or more vertices")
    end
    new(indices)
  end
end

"""
    pltype(connectivity)

Return the polytope type that the connectivity represents.
"""
pltype(::Type{Connectivity{PL,N}}) where {PL,N} = PL
pltype(c::Connectivity) = pltype(typeof(c))

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
indices(c::Connectivity) = c.indices

"""
    connect(indices, [PL])

Connect a list of `indices` from a global vector of [`Point`](@ref)
into a [`Polytope`](@ref) of type `PL`.

The type `PL` can be a [`Ngon`](@ref) in which case the length of
the indices is used to identify the actual polytope type.

Finally, the type `PL` can be ommitted. In this case, the indices are
assumed to be connected as a [`Ngon`](@ref) or as a [`Segment`](@ref).

## Examples

Connect indices into a Triangle:

```julia
connect((1,2,3), Triangle)
```

Connect indices into N-gons, a `Triangle` and a `Quadrangle`:

```julia
connect.([(1,2,3), (2,3,4,5)], Ngon)
```

Connect indices into N-gon or segment:

```julia
connect((1,2)) # Segment
connect((1,2,3)) # Triangle
connect((1,2,3,4)) # Quadrangle
```
"""
connect(indices::Tuple, PL::Type{<:Polytope}) = Connectivity{PL,length(indices)}(indices)

function connect(indices::Tuple, ::Type{Ngon})
  N = length(indices)
  Connectivity{Ngon{N},N}(indices)
end

function connect(indices::Tuple)
  N = length(indices)
  N > 2 ? connect(indices, Ngon) : connect(indices, Segment)
end

"""
    materialize(connec, points)

Materialize a face using the `connec` list and a global vector of `points`.
"""
materialize(connec::Connectivity{PL,N}, points::AbstractVector{P}) where {PL<:Polytope,N,P<:Point} =
  PL(ntuple(i -> @inbounds(points[connec.indices[i]]), N))

function Base.show(io::IO, c::Connectivity{PL}) where {PL}
  name = prettyname(PL)
  inds = c.indices
  print(io, "$name$inds")
end
