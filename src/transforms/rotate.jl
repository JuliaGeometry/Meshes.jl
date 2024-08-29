# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Rotate(R)

Rotate geometry or domain with rotation `R` from Rotations.jl.

    Rotate(u, v)

Rotation mapping the axis directed by `u` to the axis directed by `v`. 
More precisely, it maps the plane passing through the origin with normal 
vector `u` to the plane passing through the origin with normal vector `v`.

    Rotate(θ)

Rotate the 2D geometry or domain by angle `θ`, in radians, using the
`Angle2d` rotation.

## Examples

```julia
Rotate(one(RotXYZ{Float64})) # identity rotation
Rotate(AngleAxis(0.2, 1.0, 0.0, 0.0)) # rotate 0.2 radians around X axis
Rotate(rand(QuatRotation{Float64})) # random rotation
Rotate(Vec(1, 0, 0), Vec(1, 1, 1)) # rotation from (1, 0, 0) to (1, 1, 1)
Rotate(π / 2) # 2D rotation with angle in radians
```
"""
struct Rotate{R<:Rotation} <: CoordinateTransform
  rot::R
end

Rotate(u::Vec, v::Vec) = Rotate(urotbetween(u, v))

Rotate(u::Tuple, v::Tuple) = Rotate(Vec(u), Vec(v))

Rotate(θ) = Rotate(Angle2d(θ))

parameters(t::Rotate) = (; rot=t.rot)

isaffine(::Type{<:Rotate}) = true

isrevertible(::Type{<:Rotate}) = true

isinvertible(::Type{<:Rotate}) = true

inverse(t::Rotate) = Rotate(inv(t.rot))

applycoord(t::Rotate, p::Point) = withcrs(p, applycoord(t, to(p)))

applycoord(t::Rotate, v::Vec) = urotapply(t.rot, v)

# --------------
# SPECIAL CASES
# --------------

applycoord(t::Rotate, b::Box) = TransformedGeometry(b, t)

applycoord(t::Rotate, e::Ellipsoid) = Ellipsoid(radii(e), applycoord(t, center(e)), t.rot * rotation(e))

applycoord(t::Rotate, g::CartesianGrid) = TransformedGrid(g, t)

applycoord(t::Rotate, g::RectilinearGrid) = TransformedGrid(g, t)

applycoord(t::Rotate, g::StructuredGrid) = TransformedGrid(g, t)
