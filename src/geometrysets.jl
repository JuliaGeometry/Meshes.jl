# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GeometrySet(geometries)

A set of `geometries` seen as a single entity.

In a geographic map, countries can be described with
multiple polygons (a.k.a. MultiPolygon).
"""
struct GeometrySet{Dim,T,G<:Geometry{Dim,T}} <: Domain{Dim,T}
  geometries::Vector{G}
end

# -----------------
# DOMAIN INTERFACE
# -----------------

Base.getindex(gset::GeometrySet, ind::Int) = getindex(gset.geometries, ind)

nelements(gset::GeometrySet) = length(gset.geometries)
