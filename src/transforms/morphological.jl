# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Morphological(func)

Morphological transform given by a function `func`.

# Examples

```julia
TODO
f(c) = Cartesian(cospi(2*c.x), sinpi(2*c.x))
box = Box((0.0,), (1.0,))
t = Morphological(f)
tr = TransformedGeometry(box, t)

box = Box((-pi,), (pi,))
f(theta) = 1 - 1/3 * sin(2*theta)^2
tr = TransformedGeometry(box, Morphological(c -> Polar(f(c.x), c.x)))
```
"""
struct Morphological{F<:Function} <: CoordinateTransform
  func::F
  function Morphological(func::F) where {F}
    new{typeof(func)}(func)
  end
end

applycoord(t::Morphological, p::Point) = Point(t.func(coords(p)))
applycoord(t::Morphological, v::Vec) = Point(t.func(Point(v...)))
