# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    pointify(object)

Convert `object` into a vector of [`Point`](@ref).
"""
function pointify end

pointify(p::Point) = [p]

pointify(s::Sphere) = pointify(discretize(s))

pointify(p::Primitive) = pointify(boundary(p))

pointify(p::Polytope) = vertices(p)

pointify(p::PointSet) = collect(p)

pointify(m::Mesh) = vertices(m)

pointify(d::Domain) = mapreduce(pointify, vcat, d)

pointify(d::Data) = pointify(domain(d))