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
  rot = Makie.@lift TB.parameters($trans).rot
  θ = Makie.@lift first(Rotations.params($rot))

  # plot and rotate the grid
  viz!(plot, grid; color, alpha, colorscheme, segmentsize, showfacets, facetcolor)
  Makie.rotate!(plot, θ[])
end

function Makie.data_limits(plot::Viz{<:Tuple{RotatedGrid{2}}})
  tgrid = plot[:object][]
  grid = parent(tgrid)
  trans = Meshes.transform(tgrid)
  bbox = _bbox(grid, trans)
  pmin = Makie.Point3f(coordinates(minimum(bbox))..., 0)
  pmax = Makie.Point3f(coordinates(maximum(bbox))..., 0)
  Makie.limits_from_transformed_points([pmin, pmax])
end

function _bbox(grid, trans)
  q = convert(Quadrangle, boundingbox(grid))
  boundingbox(trans(q))
end
