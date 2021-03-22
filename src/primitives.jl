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

include("primitives/line.jl")
include("primitives/ray.jl")
include("primitives/plane.jl")
include("primitives/box.jl")
include("primitives/ball.jl")
include("primitives/sphere.jl")
include("primitives/cylinder.jl")
include("primitives/bezier.jl")
