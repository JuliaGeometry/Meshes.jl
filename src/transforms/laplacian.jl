# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LaplacianSmoothing(n; λ=0.5)

Perform `n` iterations of Laplacian smoothing with parameter `λ`.

## References

* Sorkine, O. 2005. [Laplacian Mesh Processing]
  (http://dx.doi.org/10.2312/egst.20051044)
"""
struct LaplacianSmoothing{T} <: StatelessGeometricTransform
  n::Int
  λ::T

  function LaplacianSmoothing{T}(n, λ) where {T}
    new(n, λ)
  end
end

function LaplacianSmoothing(n; λ=0.5)
  LaplacianSmoothing{typeof(λ)}(n, λ)
end

isrevertible(::Type{<:LaplacianSmoothing}) = true

preprocess(::LaplacianSmoothing, mesh) =
  laplacematrix(mesh, weights=:uniform)

function applypoint(transform::LaplacianSmoothing, points, prep)
  n = transform.n
  λ = transform.λ
  L = prep
  _laplacian(points, L, n, λ), L
end

function revertpoint(transform::LaplacianSmoothing, newpoints, pcache)
  n = transform.n
  λ = transform.λ
  L = pcache
  _laplacian(newpoints, L, n, λ, revert=true)
end

function _laplacian(points, L, n, λ; revert=false)
  # matrix with coordinates (nvertices x ndims)
  X = reduce(hcat, coordinates.(points)) |> transpose

  # choose between apply and revert mode
  λ₁ = revert ? -λ : λ

  # Laplace updates
  for _ in 1:n
    X = X + λ₁*L*X
  end

  # new points
  Point.(eachrow(X))
end