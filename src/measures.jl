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

function measure(c::Cylinder)
  if isright(c)
    _measure_isright(c)
  elseif hasintersectingplanes(c)
    error("Unable to calculate measure of cylinders with intersecting planes.")
  else
    c₁, c₂ = _partitioncylinder(c)
    meas(c₀) = isright(c₀) ? _measure_isright(c₀) : _measure_truncated(c₀)
    meas(c₁) + meas(c₂)
  end
end

# Get measure(c::Cylinder) for a right-cylinder (c.top ∥ c.bot)
function _measure_isright(c::Cylinder{T}) where {T}
  @assert isright(c) "This function assumes `c` is a right-cylinder."

  r = radius(c)
  ax = axis(c)
  h̄ = ax(T(1)) - ax(T(0)) # vector pointing from bottom to top of cyl axis
  h = norm(h̄)
  h * T(π) * r^2
end

# Partition a cylinder with non-parallel planes into two truncated cylinder segments
#   whose bottom planes are orthogonal to their axes
function _partitioncylinder(c::Cylinder{T}) where {T}
  @assert !isright(c) "This function assumes `c` is not a right-cylinder."

  t = top(c)
  b = bottom(c)
  r = radius(c)
  ax = axis(c)
  b̄ = ax(T(0))
  t̄ = ax(T(1))
  h̄ = t̄ - b̄ # vector pointing from bottom to top of cyl axis
  â = normalize(h̄) # unit vector pointing from bottom to top of cyl axis

  # Working from a non-right side, find a safe location for the bisecting plane
  ∠b = ∠(normal(b), â)
  if !isapprox(∠b, T(0))
    Δh = r * sin(∠b)
    bpoint = b̄ - Δh * â
  else
    ∠t = ∠(normal(t), â)
    Δh = r * sin(∠t)
    bpoint = t̄ + Δh * â
  end
  bplane = Plane(bpoint, â)

  # Construct two new cylinders by bisecting the original along bplane
  c₁ = Cylinder(bplane, b, r)
  c₂ = Cylinder(bplane, t, r)
  (c₁, c₂)
end

# Get measure(c::Cylinder) for a truncated cylinder (c.bot ⟂ axis(c))
# Assumptions: only the bottom plane is orthogonal to axis
function _measure_truncated(c::Cylinder{T}) where {T}
  r = radius(c)
  b = bottom(c)
  t = top(c)
  n̂ₜ = normal(t)
  ax = axis(c)
  â = normalize(ax(T(1)) - ax(T(0)))  # unit vector pointing from bottom to top of cyl axis

  @assert !isright(c) "This function assumes `c` is not a right-cylinder."
  @assert isapprox(normal(c.bot) × â, T(0)) "`c`s bottom plane must be orthogonal to its axis."

  # Find the Points associated with the minimum and maximum edge lengths
  t̄ = t(T(0), T(0))
  r̂ = n̂ₜ × (â × n̂ₜ)
  h̄₁ = t̄ + r * r̂
  h̄₂ = t̄ - r * r̂

  # Find edge lengths by projecting a Ray from these Points down onto the bottom plane
  h₁ = norm(h̄₁ - intersect(b, Ray(h̄₁, -â)))
  h₂ = norm(h̄₂ - intersect(b, Ray(h̄₂, -â)))

  h = T(1/2) * (h₁ + h₂)
  h * T(π) * r^2
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
