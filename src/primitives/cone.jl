# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Cone(base, apex)

A cone with `base` disk and `apex`.
See <https://en.wikipedia.org/wiki/Cone>.

See also [`ConeSurface`](@ref).
"""
struct Cone{CRS,D<:Disk{CRS},P<:Point{3,CRS}} <: Primitive{3,CRS}
  base::D
  apex::P
end

Cone(base::Disk{CRS}, apex::Tuple) where {CRS<:Cartesian} = Cone(base, Point(CRS(apex)))

paramdim(::Type{<:Cone}) = 3

base(c::Cone) = c.base

apex(c::Cone) = c.apex

height(c::Cone) = norm(center(base(c)) - apex(c))

halfangle(c::Cone) = atan(radius(base(c)), height(c))

Random.rand(rng::Random.AbstractRNG, ::Type{Cone}) = Cone(rand(rng, Disk), rand(rng, Point{3}))
