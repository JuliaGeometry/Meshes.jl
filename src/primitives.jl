# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Primitive{Dim,T}

We say that a geometry is a primitive when it can be expressed as a single
entity with no parts (a.k.a. atomic). For example, a sphere is a primitive
described in terms of a mathematical expression involving a metric and a radius.
See https://en.wikipedia.org/wiki/Geometric_primitive.
"""
abstract type Primitive{Dim,T} <: Geometry{Dim,T} end

function Base.show(io::IO, geom::G) where {G<:Primitive}
  ioctx = IOContext(io, :compact => true)
  name = prettyname(G)
  print(ioctx, "$name(")
  vals = map(fieldnames(G)) do field
    val = getfield(geom, field)
    str = repr(val, context=ioctx)
    "$field = $str"
  end
  join(ioctx, vals, ", ")
  print(ioctx, ")")
end

function Base.show(io::IO, ::MIME"text/plain", geom::G) where {G<:Primitive}
  name = prettyname(G)
  Dim = embeddim(geom)
  T = coordtype(geom)
  print(io, "$name{$Dim,$T}")
  fnames = fieldnames(G)
  len = length(fnames)
  for (i, field) in enumerate(fnames)
    div = i == len ? "\n└─ " : "\n├─ "
    val = getfield(geom, field)
    str = repr(val, context=io)
    print(io, "$div$field = $str")
  end
end

include("primitives/point.jl")
include("primitives/ray.jl")
include("primitives/line.jl")
include("primitives/bezier.jl")
include("primitives/plane.jl")
include("primitives/box.jl")
include("primitives/ball.jl")
include("primitives/sphere.jl")
include("primitives/disk.jl")
include("primitives/circle.jl")
include("primitives/cylinder.jl")
include("primitives/cylindersurface.jl")
include("primitives/cone.jl")
include("primitives/conesurface.jl")
include("primitives/torus.jl")
