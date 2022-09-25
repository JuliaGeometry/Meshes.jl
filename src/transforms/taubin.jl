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

TAPI.isrevertible(::Type{<:TaubinSmoothing}) = true

function TAPI.apply(transform::TaubinSmoothing, mesh)
  n = transform.n
  λ = transform.λ
  μ = transform.μ

  # Laplacian matrix with uniform weights
  L = laplacematrix(mesh, weights=:uniform)

  # perform Taubin updates
  newpoints = _taubin(vertices(mesh), L, n, λ, μ)

  # new mesh with same topology
  newmesh = SimpleMesh(newpoints, topology(mesh))

  newmesh, L
end

function TAPI.revert(transform::TaubinSmoothing, newmesh, cache)
  n = transform.n
  λ = transform.λ
  μ = transform.μ
  L = cache

  # reverse Taubin updates
  points = _taubin(vertices(newmesh), L, n, λ, μ, revert=true)

  # mesh with same topology
  SimpleMesh(points, topology(newmesh))
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