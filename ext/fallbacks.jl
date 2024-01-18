# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

Makie.plottype(::Geometry) = Viz{<:Tuple{Geometry}}
Makie.plottype(::Domain) = Viz{<:Tuple{Domain}}

Makie.convert_arguments(::Type{<:Viz}, geom::Geometry) = (GeometrySet([geom]),)
Makie.convert_arguments(::Type{<:Viz}, domain::Domain) = (GeometrySet(collect(domain)),)
Makie.convert_arguments(::Type{<:Viz}, mesh::Mesh) = (convert(SimpleMesh, mesh),)

# skip conversion for these types
Makie.convert_arguments(::Type{<:Viz}, gset::GeometrySet) = (gset,)
Makie.convert_arguments(::Type{<:Viz}, mesh::SimpleMesh) = (mesh,)
Makie.convert_arguments(::Type{<:Viz}, grid::Grid) = (grid,)
Makie.convert_arguments(::Type{<:Viz}, grid::SubCartesianGrid) = (grid,)
Makie.convert_arguments(::Type{<:Viz}, poly::Polygon{2}) = (poly,)

# vector of geometries for convenience
Makie.plottype(::AbstractVector{<:Geometry}) = Viz{<:Tuple{AbstractVector{<:Geometry}}}
Makie.convert_arguments(::Type{<:Viz}, geoms::AbstractVector{<:Geometry}) = (GeometrySet(geoms),)
