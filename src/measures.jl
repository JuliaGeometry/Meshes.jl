# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# ---------
# MEASURES
# ---------

"""
    measure(object)

Return the measure or "volume" of the `object`.

### Notes

- Type aliases are [`length`](@ref), [`area`](@ref), [`volume`](@ref)
"""
function measure end

measure(p::Point) = zero(coordtype(p))

measure(r::Ray) = typemax(coordtype(r))

measure(l::Line) = typemax(coordtype(l))

measure(p::Plane) = typemax(coordtype(p))

measure(b::Box) = prod(maximum(b) - minimum(b))

# https://en.wikipedia.org/wiki/Volume_of_an_n-ball
function measure(b::Ball{Dim}) where {Dim}
  r, n = radius(b), Dim
  (π^(n / 2) * r^n) / gamma(n / 2 + 1)
end

# https://en.wikipedia.org/wiki/N-sphere#Volume_and_surface_area
function measure(s::Sphere{Dim}) where {Dim}
  r, n = radius(s), Dim
  2π^(n / 2) * r^(n - 1) / gamma(n / 2)
end

function measure(d::Disk)
  T = numtype(coordtype(d))
  T(π) * radius(d)^2
end

function measure(c::Circle)
  T = numtype(coordtype(c))
  2 * T(π) * radius(c)
end

function measure(c::Cylinder)
  T = numtype(coordtype(c))
  t = top(c)
  b = bottom(c)
  r = radius(c)
  norm(t(0, 0) - b(0, 0)) * T(π) * r^2
end

function measure(c::CylinderSurface)
  T = numtype(coordtype(c))
  t = top(c)
  b = bottom(c)
  r = radius(c)
  (norm(t(0, 0) - b(0, 0)) + r) * 2 * r * T(π)
end

function measure(p::ParaboloidSurface)
  T = numtype(coordtype(c))
  r = radius(p)
  f = focallength(p)
  T(8π / 3) * f^2 * ((T(1) + r^2 / (2f)^2)^(3 / 2) - T(1))
end

# https://en.wikipedia.org/wiki/Torus
function measure(t::Torus)
  T = numtype(coordtype(c))
  R, r = radii(t)
  4T(π)^2 * R * r
end

measure(s::Segment) = norm(maximum(s) - minimum(s))

measure(t::Triangle{2}) = abs(signarea(t))

function measure(t::Triangle{3})
  A, B, C = vertices(t)
  norm((B - A) × (C - A)) / 2
end

function measure(t::Tetrahedron)
  A, B, C, D = vertices(t)
  abs((A - D) ⋅ ((B - D) × (C - D))) / 6
end

measure(c::Chain) = sum(measure, segments(c))

measure(g::Geometry) = sum(measure, simplexify(g))

measure(m::Multi) = sum(measure, parent(m))

measure(d::PointSet) = zero(coordtype(d))

measure(d::Domain) = sum(measure, d)

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

perimeter(l::Line) = zero(coordtype(l))

perimeter(p::Plane) = zero(coordtype(p))

perimeter(s::Sphere) = zero(coordtype(s))

perimeter(e::Ellipsoid) = zero(coordtype(e))
