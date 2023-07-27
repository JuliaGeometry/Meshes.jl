# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

Makie.plottype(::AbstractVector{<:Geometry}) = Viz

Makie.convert_arguments(P::Type{<:Viz}, geoms::AbstractVector{<:Geometry}) =
  Makie.convert_arguments(P, GeometrySet(geoms))

Makie.plottype(::Geometry) = Viz

Makie.convert_arguments(P::Type{<:Viz}, geom::Geometry) =
  Makie.convert_arguments(P, GeometrySet([geom]))

Makie.plottype(::Domain) = Viz

function Makie.plot!(plot::Viz{<:Tuple{Domain}})
  # retrieve domain object
  domain = plot[:object]

  # fallback to vector recipe
  viz!(plot, (Makie.@lift collect($domain)),
    color       = plot[:color],
    alpha       = plot[:alpha],
    colorscheme = plot[:colorscheme],
    facetcolor  = plot[:facetcolor],
    showfacets  = plot[:showfacets],
    pointsize   = plot[:pointsize],
    segmentsize = plot[:segmentsize]
  )
end
