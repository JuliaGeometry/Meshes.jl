# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const RotatedGrid{Dim,T} = TransformedGrid{Dim,T,G,TR} where {G,TR<:Rotate{<:Angle2d}}

function vizgrid2D!(plot::Viz{<:Tuple{RotatedGrid}})
  tgrid = plot[:object]
  grid = Makie.@lift parent(tgrid)
  trans = Makie.@lift Meshes.transform(tgrid)
  # TODO
end
