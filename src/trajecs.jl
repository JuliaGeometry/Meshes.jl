# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CylindricalTrajectory(centroids, radius)

Trajectory of cylinders of given `radius` positioned at the `centroids`.
"""
struct CylindricalTrajectory{C<:CRS,ℒ<:Len} <: Domain{C}
  centroids::Vector{Point{C}}
  radius::ℒ
  CylindricalTrajectory(centroids::Vector{Point{C}}, radius::ℒ) where {C<:CRS,ℒ<:Len} =
    new{C,float(ℒ)}(centroids, radius)
end

CylindricalTrajectory(centroids, radius::Len) = CylindricalTrajectory(collect(centroids), radius)

CylindricalTrajectory(centroids, radius) = CylindricalTrajectory(centroids, addunit(radius, u"m"))

CylindricalTrajectory(centroids::Vector{P}) where {P<:Point} = CylindricalTrajectory(centroids, oneunit(lentype(P)))

CylindricalTrajectory(centroids) = CylindricalTrajectory(collect(centroids))

topology(t::CylindricalTrajectory) = GridTopology(length(t.centroids))

function element(t::CylindricalTrajectory, ind::Int)
  ℒ = lentype(t)
  T = numtype(ℒ)
  u = unit(ℒ)
  c = t.centroids
  r = t.radius
  n = length(c)

  if n == 1 # single vertical cylinder
    p₁ = c[1] - Vec(T(0) * u, T(0) * u, T(0.5) * u)
    p₂ = c[1] + Vec(T(0) * u, T(0) * u, T(0.5) * u)
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
