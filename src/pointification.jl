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

pointify(p::Polytope) = _asvector(vertices(p))

pointify(m::Multi) = pointify(collect(m))

pointify(p::PointSet) = collect(p)

pointify(m::Mesh) = vertices(m)

pointify(d::Data) = pointify(domain(d))

# fallback with iterator of geometries
pointify(geoms) = mapreduce(pointify, vcat, geoms)

# utils
_asvector(itr) = collect(itr)
_asvector(vector::AbstractVector) = vector
