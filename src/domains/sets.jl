# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GeometrySet(geometries)

A set of `geometries` representing a [`Domain`](@ref).

## Examples

Set containing two balls centered at `(0.0, 0.0)` and `(1.0, 1.0)`:

```julia
julia> GeometrySet([Ball((0.0, 0.0)), Ball((1.0, 1.0))])
```
"""
struct GeometrySet{M<:Manifold,C<:CRS,G<:Geometry{M,C}} <: Domain{M,C}
  geoms::Vector{G}
end

# constructor with iterator of geometries
GeometrySet(geoms) = GeometrySet(map(identity, geoms))

element(d::GeometrySet, ind::Int) = d.geoms[ind]

nelements(d::GeometrySet) = length(d.geoms)

Base.parent(d::GeometrySet) = d.geoms

# specialized for efficiency
Base.vcat(d1::GeometrySet, d2::GeometrySet) = GeometrySet(vcat(d1.geoms, d2.geoms))
Base.vcat(d1::GeometrySet, d2::Domain) = GeometrySet(vcat(d1.geoms, collect(d2)))
Base.vcat(d1::Domain, d2::GeometrySet) = GeometrySet(vcat(collect(d1), d2.geoms))

# ------------------------
# SPECIAL CASE: POINT SET
# ------------------------

const PointSet{M<:Manifold,C<:CRS} = GeometrySet{M,C,Point{M,C}}

"""
    PointSet(points)

A set of `points` (a.k.a. point cloud) representing a [`Domain`](@ref).

## Examples

All point sets below are the same and contain two points in R³:

```julia
julia> PointSet([Point(1,2,3), Point(4,5,6)])
julia> PointSet(Point(1,2,3), Point(4,5,6))
julia> PointSet([(1,2,3), (4,5,6)])
julia> PointSet((1,2,3), (4,5,6))
```
"""
PointSet(points::AbstractVector{Point{M,C}}) where {M<:Manifold,C<:CRS} = PointSet{M,C}(points)
PointSet(points::Vararg{P}) where {P<:Point} = PointSet(collect(points))
PointSet(coords::AbstractVector{TP}) where {TP<:Tuple} = PointSet(Point.(coords))
PointSet(coords::Vararg{TP}) where {TP<:Tuple} = PointSet(collect(coords))

# constructor with iterator of points
PointSet(points) = PointSet(map(identity, points))
