# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Rotate(rot)

Rotate geometry or mesh with rotation `rot`
from Rotations.jl.

## Examples

```julia
Rotate(one(RotXYZ{Float64}))  # Generate identity rotation
Rotate(AngleAxis(0.2, 1.0, 0.0, 0.0))  # Rotate 0.2 radians around X-axis
Rotate(rand(QuatRotation{Float64}))  # Generate random rotation
```
"""
struct Rotate{R<:Rotation} <: CoordinateTransform
  rot::R
end

"""
    Rotate(u, v)

Rotation mapping the axis directed by `u` to the axis directed by `v`. 
More precisely, it maps the plane passing through the origin with normal 
vector `u` to the plane passing through the origin with normal vector `v`.

## Examples

```julia
Rotate(Vec(1, 0, 0), Vec(1, 1, 1))
```
"""
Rotate(u::Vec, v::Vec) = Rotate(urotbetween(u, v))

Rotate(u::Tuple, v::Tuple) = Rotate(Vec(u), Vec(v))

"""
    Rotate(θ)

Rotate the 2D geometry or mesh by angle `θ`, in radians,
using the `Angle2d` rotation.

## Examples

```julia
Rotate(π / 2)
```
"""
Rotate(θ) = Rotate(Angle2d(θ))

parameters(t::Rotate) = (; rot=t.rot)

isaffine(::Type{<:Rotate}) = true

isrevertible(::Type{<:Rotate}) = true

isinvertible(::Type{<:Rotate}) = true

inverse(t::Rotate) = Rotate(inv(t.rot))

applycoord(t::Rotate, v::Vec) = urotapply(t.rot, v)

# --------------
# SPECIAL CASES
# --------------

applycoord(t::Rotate, b::Box) = _applycoord(t, b, Val(embeddim(b)))

_applycoord(t::Rotate, b::Box, ::Val{2}) = applycoord(t, convert(Quadrangle, b))

_applycoord(t::Rotate, b::Box, ::Val{3}) = applycoord(t, convert(Hexahedron, b))

applycoord(t::Rotate, g::CartesianGrid) = TransformedGrid(g, t)

applycoord(t::Rotate, g::RectilinearGrid) = TransformedGrid(g, t)

applycoord(t::Rotate, g::StructuredGrid) = TransformedGrid(g, t)
