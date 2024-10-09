# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Morphological(func)

Morphological transform given by a function `func`.

# Examples

```julia
TODO
ball = Ball((0.0, 0.0), 1)
tr = TransformedGeometry(ball, Morphological(c -> Cartesian(c.x + c.y, c.y, c.x - c.y)))
```
"""
struct Morphological{F<:Function} <: CoordinateTransform
  func::F
  function Morphological(func::F) where {F}
    new{typeof(func)}(func)
  end
end

applycoord(t::Morphological, p::Point) = Point(t.func(coords(p)))
