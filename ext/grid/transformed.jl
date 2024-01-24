# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const RotatedGrid{Dim,T} = TransformedGrid{Dim,T,G,TR} where {G,TR<:Rotate{<:Angle2d}}

function vizgrid2D!(plot::Viz{<:Tuple{RotatedGrid}})
  tgrid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]
  segmentsize = plot[:segmentsize]
  showfacets = plot[:showfacets]
  facetcolor = plot[:facetcolor]

  grid = Makie.@lift parent($tgrid)
  trans = Makie.@lift Meshes.transform($tgrid)
  rot = Makie.@lift first(TB.parameters($trans))
  θ = Makie.@lift first(Rotations.params($rot))

  # plot and rotate the grid
  viz!(plot, grid; color, alpha, colorscheme, segmentsize, showfacets, facetcolor)
  Makie.rotate!(plot, θ[])
end

const TranslatedGrid{Dim,T} = TransformedGrid{Dim,T,G,TR} where {G,TR<:Translate{Dim}}

vizgrid2D!(plot::Viz{<:Tuple{TranslatedGrid}}) = translatedgrid!(plot)

vizgrid3D!(plot::Viz{<:Tuple{TranslatedGrid}}) = translatedgrid!(plot)

function translatedgrid!(plot)
  tgrid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]
  segmentsize = plot[:segmentsize]
  showfacets = plot[:showfacets]
  facetcolor = plot[:facetcolor]

  grid = Makie.@lift parent($tgrid)
  trans = Makie.@lift Meshes.transform($tgrid)
  offsets = Makie.@lift first(TB.parameters($trans))

  # plot and translate the grid
  viz!(plot, grid; color, alpha, colorscheme, segmentsize, showfacets, facetcolor)
  Makie.translate!(plot, offsets[]...)
end
