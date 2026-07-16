# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# ---------
# MEASURES
# ---------

"""
    measure(object)

Return the measure of the geometric `object`.

This function is also known as [`length`](@ref),
[`area`](@ref) or [`volume`](@ref) depending on
the parametric dimension of the object.
"""
measure(g::GeometryOrDomain) = integral(_ -> 1, g)

measure(p::Point) = zero(lentype(p))

measure(r::Ray) = typemax(lentype(r))

measure(l::Line) = typemax(lentype(l))

measure(p::Plane) = typemax(lentype(p))^2

measure(b::Box{<:𝔼}) = prod(maximum(b) - minimum(b))

# https://en.wikipedia.org/wiki/Volume_of_an_n-ball
function measure(b::Ball{<:𝔼})
  T = numtype(lentype(b))
  r, n = radius(b), embeddim(b)
  T(π)^T(n / 2) * r^n / gamma(T(n / 2) + 1)
end

# https://en.wikipedia.org/wiki/N-sphere#Volume_and_surface_area
function measure(s::Sphere{<:𝔼})
  T = numtype(lentype(s))
  r, n = radius(s), embeddim(s)
  2 * T(π)^T(n / 2) * r^(n - 1) / gamma(T(n / 2))
end

measure(d::Disk) = π * radius(d)^2

measure(c::Circle) = 2 * (π * radius(c))

function measure(c::Cylinder)
  t = top(c)
  b = bottom(c)
  r = radius(c)
  h = norm(t(0, 0) - b(0, 0))
  π * r^2 * h
end

function measure(c::CylinderSurface)
  t = top(c)
  b = bottom(c)
  r = radius(c)
  h = norm(t(0, 0) - b(0, 0))
  2 * (π * r) * (h + r)
end

function measure(p::ParaboloidSurface)
  T = numtype(lentype(p))
  r = radius(p)
  f = focallength(p)
  (8 * T(π) / 3) * f^2 * ((1 + r^2 / (2f)^2)^T(3 / 2) - 1)
end

# https://en.wikipedia.org/wiki/Torus
function measure(t::Torus)
  T = numtype(lentype(t))
  R, r = radii(t)
  4 * T(π)^2 * R * r
end

measure(s::Segment{<:𝔼}) = norm(maximum(s) - minimum(s))

measure(c::Chain{<:𝔼}) = sum(measure, segments(c))

function measure(t::Triangle{<:𝔼})
  A, B, C = vertices(t)
  norm((B - A) × (C - A)) / 2
end

function measure(t::Tetrahedron)
  A, B, C, D = vertices(t)
  abs((A - D) ⋅ ((B - D) × (C - D))) / 6
end

function measure(p::Polygon{𝔼{2}})
  Σ = sum(rings(p)) do r
    v = vertices(r)
    n = nvertices(r)
    o = coordmean(v)
    sum(signarea(o, v[i], v[i + 1]) for i in 1:n)
  end
  abs(Σ)
end

measure(m::Multi{<:𝔼}) = sum(measure, parent(m))

measure(g::Geometry{<:𝔼}) = sum(measure, simplexify(g))

measure(d::PointSet) = zero(lentype(d))

measure(d::OrthoRegularGrid) = nelements(d) * measure(first(d))

measure(d::Domain{<:𝔼}) = sum(measure, d)

# --------
# ALIASES
# --------

"""
    length(object)

Return the length of the `object`.

See also [`measure`](@ref).
"""
function Base.length(g::GeometryOrDomain)
  if isparametrized(g) && paramdim(g) != 1
    throw(ArgumentError("invalid parametric dimension for computing length"))
  end
  measure(g)
end

"""
    area(object)

Return the area of the `object`.

See also [`measure`](@ref).
"""
function area(g::GeometryOrDomain)
  if isparametrized(g) && paramdim(g) != 2
    throw(ArgumentError("invalid parametric dimension for computing area"))
  end
  measure(g)
end

"""
    volume(object)

Return the volume of the `object`.

See also [`measure`](@ref).
"""
function volume(g::GeometryOrDomain)
  if isparametrized(g) && paramdim(g) != 3
    throw(ArgumentError("invalid parametric dimension for computing volume"))
  end
  measure(g)
end

# ----------
# PERIMETER
# ----------

"""
    perimeter(object)

Return the perimeter of the `object`, i.e.
the [`measure`](@ref) of its [`boundary`](@ref).
"""
perimeter(g) = measure(boundary(g))

perimeter(s::Segment) = zero(lentype(s))

perimeter(c::BezierCurve) = zero(lentype(c))

perimeter(c::ParametrizedCurve) = zero(lentype(c))

perimeter(r::Ray) = zero(lentype(r))

perimeter(l::Line) = zero(lentype(l))

perimeter(p::Plane) = zero(lentype(p))

perimeter(s::Sphere) = zero(lentype(s))

perimeter(e::Ellipsoid) = zero(lentype(e))
