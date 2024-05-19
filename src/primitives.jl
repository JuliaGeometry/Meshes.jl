# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Primitive{Dim}

We say that a geometry is a primitive when it can be expressed as a single
entity with no parts (a.k.a. atomic). For example, a sphere is a primitive
described in terms of a mathematical expression involving a metric and a radius.
See <https://en.wikipedia.org/wiki/Geometric_primitive>.
"""
abstract type Primitive{Dim} <: Geometry{Dim} end

function Base.show(io::IO, geom::Primitive)
  name = prettyname(geom)
  ioctx = IOContext(io, :compact => true)
  print(io, "$name(")
  printfields(ioctx, geom, singleline=true)
  print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", geom::Primitive)
  summary(io, geom)
  printfields(io, geom)
end

include("primitives/point.jl")
include("primitives/ray.jl")
include("primitives/line.jl")
include("primitives/bezier.jl")
include("primitives/plane.jl")
include("primitives/box.jl")
include("primitives/ball.jl")
include("primitives/sphere.jl")
include("primitives/ellipsoid.jl")
include("primitives/disk.jl")
include("primitives/circle.jl")
include("primitives/cylinder.jl")
include("primitives/cylindersurface.jl")
include("primitives/cone.jl")
include("primitives/conesurface.jl")
include("primitives/frustum.jl")
include("primitives/frustumsurface.jl")
include("primitives/paraboloidsurface.jl")
include("primitives/torus.jl")
