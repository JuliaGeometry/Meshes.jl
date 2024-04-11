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

measure(::Point{Dim,T}) where {Dim,T} = zero(T)

measure(::Ray{Dim,T}) where {Dim,T} = typemax(T)

measure(::Line{Dim,T}) where {Dim,T} = typemax(T)

measure(::Plane{T}) where {T} = typemax(T)

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

measure(d::Disk{T}) where {T} = T(π) * radius(d)^2

measure(c::Circle{T}) where {T} = 2 * T(π) * radius(c)

function measure(c::Cylinder{T}) where {T}
  t = top(c)
  b = bottom(c)
  r = radius(c)
  norm(t(0, 0) - b(0, 0)) * T(π) * r^2
end

function measure(c::CylinderSurface{T}) where {T}
  t = top(c)
  b = bottom(c)
  r = radius(c)
  (norm(t(0, 0) - b(0, 0)) + r) * 2 * r * T(π)
end

function measure(p::ParaboloidSurface{T}) where {T}
  r = radius(p)
  f = focallength(p)
  T(8π / 3) * f^2 * ((T(1) + r^2 / (2f)^2)^(3 / 2) - T(1))
end

# https://en.wikipedia.org/wiki/Torus
function measure(t::Torus{T}) where {T}
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

measure(::PointSet{Dim,T}) where {Dim,T} = zero(T)

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

perimeter(::Line{Dim,T}) where {Dim,T} = zero(T)

perimeter(::Plane{T}) where {T} = zero(T)

perimeter(::Sphere{Dim,T}) where {Dim,T} = zero(T)

perimeter(::Ellipsoid{T}) where {T} = zero(T)
