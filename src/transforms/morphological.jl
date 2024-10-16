# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Morphological(fun; boundary=false)

Morphological transform given by a function `fun`
that maps the coordinates of a geometry or a domain
to new coordinates (`coords -> newcoords`).

# Examples

```julia
ball = Ball((0, 0), 1)
ball |> Morphological(c -> Cartesian(c.x + c.y, c.y, c.x - c.y))

triangle = Triangle(latlon(0, 0), latlon(0, 45), latlon(45, 0))
triangle |> Morphological(c -> LatLonAlt(c.lon, c.lat, 0.0m), boundary=true)
```
"""
struct Morphological{Boundary,F<:Function} <: CoordinateTransform
  fun::F
  Morphological{Boundary}(fun::F) where {Boundary,F<:Function} = new{Boundary,F}(fun)
end

Morphological(fun; boundary=false) = Morphological{boundary}(fun)

parameters(t::Morphological{Boundary}) where {Boundary} = (fun=t.fun, boundary=Boundary)

applycoord(t::Morphological, p::Point) = Point(t.fun(coords(p)))

applycoord(::Morphological, v::Vec) = v

# --------------
# SPECIAL CASES
# --------------

applycoord(t::Morphological, g::Primitive) = TransformedGeometry(g, t)

applycoord(t::Morphological{true}, g::Polytope) = TransformedGeometry(g, t)

applycoord(t::Morphological, g::RegularGrid) = TransformedGrid(g, t)

applycoord(t::Morphological, g::RectilinearGrid) = TransformedGrid(g, t)

applycoord(t::Morphological, g::StructuredGrid) = TransformedGrid(g, t)

# -----------
# IO METHODS
# -----------

Base.show(io::IO, t::Morphological{Boundary}) where {Boundary} =
  print(io, "Morphological(fun: $(t.fun), boundary: $Boundary)")

function Base.show(io::IO, ::MIME"text/plain", t::Morphological{Boundary}) where {Boundary}
  summary(io, t)
  println(io)
  println(io, "├─ fun: $(t.fun)")
  print(io, "└─ boundary: $Boundary")
end
