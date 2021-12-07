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
struct TaubinSmoothing <: SmoothingMethod
  n::Int
  λ::Float64
  μ::Float64
end

function TaubinSmoothing(n; λ=0.5, μ=-0.5)
  @assert n > 0 "invalid number of iterations"
  @assert 0 < λ ≤ -μ < 1 "invalid parameters λ and μ"
  TaubinSmoothing(n, λ, μ)
end

function smooth(mesh, method::TaubinSmoothing)
  n = method.n
  λ = method.λ
  μ = method.μ

  # Laplacian matrix with uniform weights
  L = laplacematrix(mesh, weights=:uniform)

  # matrix with vertex coordinates (nvertices x ndims)
  V = reduce(hcat, coordinates.(vertices(mesh))) |> transpose

  # Taubin updates
  for _ in 1:n
    V = V + λ*L*V
    V = V + μ*L*V
  end

  # new points of the smooth mesh
  points = Point.(eachrow(V))

  SimpleMesh(points, topology(mesh))
end
