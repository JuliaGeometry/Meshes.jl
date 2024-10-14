# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Morphological(fun)

Morphological transform given by a function `fun`
that maps the coordinates of a geometry or a domain
to new coordinates (`coords -> newcoords`).

# Examples

```julia
ball = Ball((0, 0), 1)
ball |> Morphological(c -> Cartesian(c.x + c.y, c.y, c.x - c.y))
```
"""
struct Morphological{F<:Function} <: CoordinateTransform
  fun::F
end

parameters(t::Morphological) = (; fun=t.fun)

applycoord(t::Morphological, p::Point) = Point(t.fun(coords(p)))

applycoord(::Morphological, v::Vec) = v

# --------------
# SPECIAL CASES
# --------------

# method to fix ambiguities
applycoord(t::Morphological, g::TransformedGeometry) = TransformedGeometry(g, t)

applycoord(t::Morphological, g::Geometry) = TransformedGeometry(g, t)

applycoord(t::Morphological, g::RegularGrid) = TransformedGrid(g, t)

applycoord(t::Morphological, g::RectilinearGrid) = TransformedGrid(g, t)

applycoord(t::Morphological, g::StructuredGrid) = TransformedGrid(g, t)
