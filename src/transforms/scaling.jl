# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ScaleCoords([s1, s2, ...])
    ScaleCoords(s1, s2, ...)

Transform geometry or mesh by scaling coordinates
`(x1, x2, ...) ↦ (s1*x1, s2*x2, ...)`. 
"""
struct ScaleCoords{T} <: GeometricTransform
  factors::Vector{T}
  function ScaleCoords{T}(factors) where {T}
    if any(≤(0), factors)
      throw(ArgumentError("Scaling factors must be positive."))
    end
    new(factors)
  end
end

ScaleCoords(factors::AbstractVector) =
  ScaleCoords{eltype(factors)}(factors)

ScaleCoords(factors...) =
  ScaleCoords(collect(factors))

isrevertible(::Type{<:ScaleCoords}) = true

preprocess(transform::ScaleCoords, object) = transform.factors

function applypoint(::ScaleCoords, points, prep)
  s = prep
  newpoints = [Point(s .* coordinates(p)) for p in points]
  newpoints, prep
end

function revertpoint(::ScaleCoords, newpoints, cache)
  s = cache
  [Point((1 ./ s) .* coordinates(p)) for p in newpoints]
end

function reapplypoint(::ScaleCoords, points, cache)
  s = cache
  [Point(s .* coordinates(p)) for p in points]
end
