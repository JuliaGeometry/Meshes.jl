# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CylindricalTrajectory(centroids, radius)

Trajectory of cylinders of given `radius` positioned at the `centroids`.
"""
struct CylindricalTrajectory{T} <: Domain{3,T}
  centroids::Vector{Point{3,T}}
  radius::T
end

CylindricalTrajectory(centroids::AbstractVector{Point{3,T}}, radius) where {T} =
  CylindricalTrajectory(centroids, T(radius))

CylindricalTrajectory(centroids) = CylindricalTrajectory(centroids, 1)

topology(t::CylindricalTrajectory) = GridTopology(length(t.centroids))

function element(t::CylindricalTrajectory{T}, ind::Int) where {T}
  c = t.centroids
  r = t.radius
  n = length(c)

  if n == 1 # single vertical cylinder
    p₁ = c[1] - Vec{3,T}(0, 0, 0.5)
    p₂ = c[1] + Vec{3,T}(0, 0, 0.5)
    return Cylinder(p₁, p₂, r)
  end

  if ind == 1 # head of trajectory
    # points at cylinder planes
    p₂ = center(Segment(c[ind], c[ind + 1]))
    p₁ = p₂ - 2 * (p₂ - c[ind])

    # normals to cylinder planes
    n₂ = c[ind + 1] - c[ind]
    n₁ = n₂
  elseif ind == n # tail of trajectory
    # points at cylinder planes
    p₁ = center(Segment(c[ind - 1], c[ind]))
    p₂ = p₁ + 2 * (c[ind] - p₁)

    # normals to cylinder planes
    n₁ = c[ind] - c[ind - 1]
    n₂ = n₁
  else # middle of trajectory
    # points at cylinder planes
    p₁ = center(Segment(c[ind - 1], c[ind]))
    p₂ = center(Segment(c[ind], c[ind + 1]))

    # normals to cylinder planes
    n₁ = c[ind] - c[ind - 1]
    n₂ = c[ind + 1] - c[ind]
  end

  # cylinder with given radius and planes
  Cylinder(Plane(p₁, n₁), Plane(p₂, n₂), r)
end

nelements(t::CylindricalTrajectory) = length(t.centroids)

Base.eltype(t::CylindricalTrajectory) = typeof(first(t))

radius(t::CylindricalTrajectory) = t.radius
