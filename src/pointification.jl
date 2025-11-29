# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    pointify(geometry)

Return vector of [`Point`](@ref)s on the
[`embedboundary`](@ref) of the `geometry`.
"""
pointify(g::Geometry) = _points(embedboundary(g))

# discretize boundary
_points(g::Geometry) = vertices(discretize(g))

# skip discretization
_points(p::Point) = [p]
_points(m::MultiPoint) = parent(m)
_points(m::Mesh) = vertices(m)
