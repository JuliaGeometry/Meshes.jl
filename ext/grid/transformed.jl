# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

isimplemented(::TB.Transform) = false
isimplemented(::Rotate{<:Angle2d}) = true
isimplemented(::Translate) = true

vizgrid2D!(plot::Viz{<:Tuple{TransformedGrid}}) = transformedgrid!(plot, vizmesh2D!)

vizgrid3D!(plot::Viz{<:Tuple{TransformedGrid}}) = transformedgrid!(plot, vizmesh3D!)

function transformedgrid!(plot, fallback)
  tgrid = plot[:object]
  grid = Makie.@lift parent($tgrid)
  trans = Makie.@lift Meshes.transform($tgrid)
  if isaffine(trans[]) && isimplemented(trans[])
    color = plot[:color]
    alpha = plot[:alpha]
    colorscheme = plot[:colorscheme]
    segmentsize = plot[:segmentsize]
    showfacets = plot[:showfacets]
    facetcolor = plot[:facetcolor]
    viz!(plot, grid; color, alpha, colorscheme, segmentsize, showfacets, facetcolor)
    makietransform!(plot, trans[])
  else
    fallback(tgrid)
  end
end

function makietransform!(plot, trans::Rotate{<:Angle2d})
  rot = first(TB.parameters(trans))
  θ = first(Rotations.params(rot))
  Makie.rotate!(plot, θ)
end

function makietransform!(plot, trans::Translate)
  offsets = first(TB.parameters(trans))
  Makie.translate!(plot, offsets...)
end
