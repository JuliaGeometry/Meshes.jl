# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LambdaMuSmoothing(n, λ, μ)

Perform `n` smoothing iterations with parameters `λ` and `μ`.

See also [`LaplaceSmoothing`](@ref), [`TaubinSmoothing`](@ref).

## References

* Taubin, G. 1995. [Curve and Surface Smoothing without Shrinkage]
  (https://ieeexplore.ieee.org/document/466848)
"""
struct LambdaMuSmoothing{T} <: StatelessGeometricTransform
  n::Int
  λ::T
  μ::T
end

isrevertible(::Type{<:LambdaMuSmoothing}) = true

preprocess(::LambdaMuSmoothing, mesh) = laplacematrix(mesh; weights=:uniform)

function applypoint(transform::LambdaMuSmoothing, points, prep)
  n = transform.n
  λ = transform.λ
  μ = transform.μ
  L = prep
  return _smooth(points, L, n, λ, μ), L
end

function revertpoint(transform::LambdaMuSmoothing, newpoints, pcache)
  n = transform.n
  λ = transform.λ
  μ = transform.μ
  L = pcache
  return _smooth(newpoints, L, n, λ, μ; revert=true)
end

function _smooth(points, L, n, λ, μ; revert=false)
  # matrix with coordinates (nvertices x ndims)
  X = transpose(reduce(hcat, coordinates.(points)))

  # choose between apply and revert mode
  λ₁, λ₂ = revert ? (-μ, -λ) : (λ, μ)

  # Taubin updates
  for _ in 1:n
    X = X + λ₁ * L * X
    X = X + λ₂ * L * X
  end

  # new points
  return Point.(eachrow(X))
end

"""
LaplaceSmoothing(n, λ=0.5)

Perform `n` iterations of Laplace smoothing with parameter `λ`.

## References

* Sorkine, O. 2005. [Laplacian Mesh Processing]
  (http://dx.doi.org/10.2312/egst.20051044)
"""
LaplaceSmoothing(n, λ=0.5) = LambdaMuSmoothing(n, λ, zero(λ))

"""
TaubinSmoothing(n, λ=0.5)

Perform `n` iterations of Taubin smoothing with parameter `0 < λ < 1`.

## References

* Taubin, G. 1995. [Curve and Surface Smoothing without Shrinkage]
  (https://ieeexplore.ieee.org/document/466848)
"""
TaubinSmoothing(n, λ=0.5) = LambdaMuSmoothing(n, λ, -λ)
