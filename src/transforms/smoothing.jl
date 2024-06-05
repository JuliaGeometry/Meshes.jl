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
struct LambdaMuSmoothing{T} <: GeometricTransform
  n::Int
  λ::T
  μ::T
end

parameters(t::LambdaMuSmoothing) = (n=t.n, λ=t.λ, μ=t.μ)

isrevertible(::Type{<:LambdaMuSmoothing}) = true

function apply(transform::LambdaMuSmoothing, mesh::Mesh)
  n = transform.n
  λ = transform.λ
  μ = transform.μ
  L = _laplacian(mesh)
  _smooth(mesh, L, n, λ, μ), L
end

function revert(transform::LambdaMuSmoothing, mesh::Mesh, cache)
  n = transform.n
  λ = transform.λ
  μ = transform.μ
  L = cache
  _smooth(mesh, L, n, λ, μ, revert=true)
end

_laplacian(mesh) = laplacematrix(mesh, kind=:uniform)

function _smooth(mesh, L, n, λ, μ; revert=false)
  # retrieve vertices
  points = vertices(mesh)

  # matrix with coordinates (nvertices x ndims)
  X = reduce(hcat, to.(points)) |> transpose

  # choose between apply and revert mode
  λ₁, λ₂ = revert ? (-μ, -λ) : (λ, μ)

  # Taubin updates
  for _ in 1:n
    X = X + λ₁ * L * X
    X = X + λ₂ * L * X
  end

  # new points
  newpoints = Point.(Tuple.(eachrow(X)))

  # new mesh
  SimpleMesh(newpoints, topology(mesh))
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
