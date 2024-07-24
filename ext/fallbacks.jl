# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

Makie.plottype(::Geometry) = Viz{<:Tuple{Geometry}}
Makie.plottype(::Domain) = Viz{<:Tuple{Domain}}
Makie.plottype(::Vec) = Viz{<:Tuple{Vec}}
Makie.plottype(::AbstractVector{<:Vec}) = Viz{<:Tuple{AbstractVector{<:Vec}}}

Makie.convert_arguments(::Type{<:Viz}, geom::Geometry) = (GeometrySet([geom]),)
Makie.convert_arguments(::Type{<:Viz}, domain::Domain) = (GeometrySet(collect(domain)),)
Makie.convert_arguments(::Type{<:Viz}, vec::Vec) = ([vec],)

# skip conversion for these types
Makie.convert_arguments(::Type{<:Viz}, mesh::Mesh) = (mesh,)
Makie.convert_arguments(::Type{<:Viz}, gset::GeometrySet) = (gset,)
Makie.convert_arguments(::Type{<:Viz}, subdom::SubDomain) = (subdom,)
Makie.convert_arguments(::Type{<:Viz}, vecs::AbstractVector{<:Vec}) = (vecs,)

# vector of geometries for convenience
Makie.plottype(::AbstractVector{<:Geometry}) = Viz{<:Tuple{AbstractVector{<:Geometry}}}
Makie.convert_arguments(::Type{<:Viz}, geoms::AbstractVector{<:Geometry}) = (GeometrySet(geoms),)
