# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

mayberepeat(value::AbstractVector, meshes) = [value[e] for (e, mesh) in enumerate(meshes) for _ in 1:nelements(mesh)]

mayberepeat(value, meshes) = value

function vizmany!(plot, meshes)
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]
  pointsize = plot[:pointsize]
  segmentsize = plot[:segmentsize]

  mesh = Makie.@lift reduce(merge, $meshes)
  colors = Makie.@lift mayberepeat($color, $meshes)
  alphas = Makie.@lift mayberepeat($alpha, $meshes)

  viz!(
    plot,
    mesh,
    color=colors,
    alpha=alphas,
    colorscheme=colorscheme,
    showfacets=false,
    pointsize=pointsize,
    segmentsize=segmentsize
  )
end
