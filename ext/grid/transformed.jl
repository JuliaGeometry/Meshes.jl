# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

isoptimized(::TB.Identity) = true
isoptimized(t::TB.SequentialTransform) = all(isoptimized, t)

isoptimized(::GeometricTransform) = false
isoptimized(::Rotate{<:Angle2d}) = true
isoptimized(::Translate) = true
isoptimized(::Scale) = true

vizgrid2D!(plot::Viz{<:Tuple{TransformedGrid}}) = transformedgrid!(plot, vizmesh2D!)

vizgrid3D!(plot::Viz{<:Tuple{TransformedGrid}}) = transformedgrid!(plot, vizmesh3D!)

function transformedgrid!(plot, fallback)
  tgrid = plot[:object]
  grid = Makie.@lift parent($tgrid)
  trans = Makie.@lift Meshes.transform($tgrid)
  if isoptimized(trans[])
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

makietransform!(plot, trans::TB.Identity) = nothing

makietransform!(plot, trans::TB.SequentialTransform) = foreach(t -> makietransform!(plot, t), trans)

function makietransform!(plot, trans::Rotate{<:Angle2d})
  rot = first(TB.parameters(trans))
  θ = first(Rotations.params(rot))
  Makie.rotate!(plot, θ)
end

function makietransform!(plot, trans::Translate)
  offsets = first(TB.parameters(trans))
  Makie.translate!(plot, offsets...)
end

function makietransform!(plot, trans::Scale)
  factors = first(TB.parameters(trans))
  Makie.scale!(plot, factors...)
end
