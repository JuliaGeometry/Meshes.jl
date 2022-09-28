# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TaubinSmoothing(n; λ=0.5, μ=-0.5)

Perform `n` iterations of Taubin smoothing with parameters `λ` and `μ`.

## References

* Taubin, G. 1995. [Curve and Surface Smoothing without Shrinkage]
  (https://ieeexplore.ieee.org/document/466848)
"""
struct TaubinSmoothing{T} <: StatelessGeometricTransform
  n::Int
  λ::T
  μ::T

  function TaubinSmoothing{T}(n, λ, μ) where {T}
    @assert n > 0 "invalid number of iterations"
    @assert 0 < λ ≤ -μ < 1 "invalid parameters λ and μ"
    new(n, λ, μ)
  end
end

function TaubinSmoothing(n; λ=0.5, μ=-0.5)
  λ′, μ′ = promote(λ, μ)
  TaubinSmoothing{typeof(λ′)}(n, λ, μ)
end

isrevertible(::Type{<:TaubinSmoothing}) = true

preprocess(::TaubinSmoothing, mesh) =
  laplacematrix(mesh, weights=:uniform)

function applypoint(transform::TaubinSmoothing, points, prep)
  n = transform.n
  λ = transform.λ
  μ = transform.μ
  L = prep
  _taubin(points, L, n, λ, μ), L
end

function revertpoint(transform::TaubinSmoothing, newpoints, pcache)
  n = transform.n
  λ = transform.λ
  μ = transform.μ
  L = pcache
  _taubin(newpoints, L, n, λ, μ, revert=true)
end

function _taubin(points, L, n, λ, μ; revert=false)
  # matrix with coordinates (nvertices x ndims)
  X = reduce(hcat, coordinates.(points)) |> transpose

  # choose between apply and revert mode
  λ₁, λ₂ = revert ? (-μ, -λ) : (λ, μ)

  # Taubin updates
  for _ in 1:n
    X = X + λ₁*L*X
    X = X + λ₂*L*X
  end

  # new points
  Point.(eachrow(X))
end