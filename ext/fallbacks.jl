# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

Makie.plottype(::AbstractVector{<:Geometry}) = Viz{<:Tuple{AbstractVector{<:Geometry}}}

Makie.convert_arguments(::Type{<:Viz}, geoms::AbstractVector{<:Geometry}) = (GeometrySet(geoms),)

Makie.plottype(::Geometry) = Viz{<:Tuple{Geometry}}

Makie.convert_arguments(::Type{<:Viz}, geom::Geometry) = (GeometrySet([geom]),)

Makie.plottype(::Domain) = Viz{<:Tuple{Domain}}

function Makie.plot!(plot::Viz{<:Tuple{Domain}})
  # retrieve domain object
  domain = plot[:object]

  # fallback to vector recipe
  viz!(
    plot,
    (Makie.@lift collect($domain)),
    color=plot[:color],
    alpha=plot[:alpha],
    colorscheme=plot[:colorscheme],
    pointsize=plot[:pointsize],
    segmentsize=plot[:segmentsize],
    showfacets=plot[:showfacets],
    facetcolor=plot[:facetcolor]
  )
end
