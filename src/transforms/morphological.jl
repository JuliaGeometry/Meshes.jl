# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Morphological(func)

Morphological transform given by a function `func`.

# Examples

```julia
TODO
# with `Cartesian`:
# f(c) = Point(cospi(2*Meshes.ustrip(c.x)), sinpi(2*Meshes.ustrip(c.x)))
# with `Vec`:
f(v) = cospi(2*v[1]), sinpi(2*v[1])
box = Box((0.0,), (1.0,))
t = Morphological(f)
tr = TransformedGeometry(box, t)

box = Box((-pi,), (pi,))
f(theta) = 1 - 1/3 * sin(2*theta)^2
# with `Cartesian`:
# tr = TransformedGeometry(box, Morphological(c -> Point(Polar(f(Meshes.ustrip(c.x)), Meshes.ustrip(c.x)))))
# with `Vec:
tr = TransformedGeometry(box, Morphological(v -> Polar(f(v[1]), v[1])))
```
"""
struct Morphological{F<:Function} <: CoordinateTransform
  func::F
  function Morphological(func::F) where {F}
    new{typeof(func)}(func)
  end
end

# applycoord(t::Morphological, p::Point) = t.func(coords(p))
# applycoord(t::Morphological, v::Vec) = t.func(Point(v...))
applycoord(t::Morphological, p::Point) = Point(t.func(ustrip(to(p))))
applycoord(t::Morphological, v::Vec) = Point(t.func(v))
