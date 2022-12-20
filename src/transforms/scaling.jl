# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Scaling(sx, sy, sz)

Transform geometry or mesh by scaling each coordinate of the vertices by 
the corresponding factor. 

## Arguments

- `sx`: scaling factor for the first coordinate
- `sy`: scaling factor for the second coordinate
- `sz`: scaling factor for the third coordinate

## Examples

```julia
Scaling(1, 2, 3)
```
"""
struct Scaling{T<:Real} <: GeometricTransform
  sx::T
  sy::T
  sz::T
end

function Scaling(sx::T1, sy::T2, sz::T3) where {T1<:Real,T2<:Real,T3<:Real}
  if sx <= 0 || sy <= 0 || sz <= 0
    error("The scaling factors must be nonnegative.")
  end
  T = promote(sx, sy, sz)
  Scaling(T(sx), T(sy), T(sz))
end

isrevertible(::Type{<:Scaling}) = true

function preprocess(transform::Scaling, object)
  [transform.sx, transform.sy, transform.sz]
end

function applypoint(::Scaling, points, prep)
  scale_vector = prep
  newpoints = [Point(scale_vector .* coordinates(p)) for p in points]
  newpoints, prep
end

function revertpoint(::Scaling, newpoints, cache)
  scale_vector = cache
  iscale_vector = 1 ./ scale_vector
  [Point(iscale_vector .* coordinates(p)) for p in newpoints]
end

function reapplypoint(::Scaling, points, cache)
  scale_vector = cache
  [Point(scale_vector .* coordinates(p)) for p in points]
end
