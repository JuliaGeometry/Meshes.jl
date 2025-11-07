# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# generic fallback for geometry
Makie.plottype(::Geometry) = Viz{<:Tuple{Geometry}}
Makie.convert_arguments(::Type{<:Viz}, geom::Geometry) = (GeometrySet([geom]),)

# generic fallback for domain
Makie.plottype(::Domain) = Viz{<:Tuple{Domain}}
Makie.convert_arguments(::Type{<:Viz}, domain::Domain) = (GeometrySet(collect(domain)),)

# skip conversion and use specialized methods
Makie.convert_arguments(::Type{<:Viz}, mesh::Mesh) = (mesh,)
Makie.convert_arguments(::Type{<:Viz}, gset::GeometrySet) = (gset,)
Makie.convert_arguments(::Type{<:Viz}, subdom::SubDomain) = (subdom,)

# vector of geometries for convenience
Makie.plottype(::AbstractVector{<:Geometry}) = Viz{<:Tuple{AbstractVector{<:Geometry}}}
Makie.convert_arguments(::Type{<:Viz}, geoms::AbstractVector{<:Geometry}) = (GeometrySet(geoms),)

# geometric vectors for convenience
Makie.plottype(::Vec) = Viz{<:Tuple{Vec}}
Makie.convert_arguments(::Type{<:Viz}, vec::Vec) = (GeometrySet([asray(vec)]),)
Makie.plottype(::AbstractVector{<:Vec}) = Viz{<:Tuple{AbstractVector{<:Vec}}}
Makie.convert_arguments(::Type{<:Viz}, vecs::AbstractVector{<:Vec}) = (GeometrySet(asray.(vecs)),)
