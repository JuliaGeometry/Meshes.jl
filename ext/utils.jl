# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

mayberepeat(value, objs) = value

mayberepeat(value::AbstractVector, objs) = [value[e] for (e, obj) in enumerate(objs) for _ in 1:length(obj)]

concat(obj₁, obj₂) = vcat(obj₁, obj₂)

concat(mesh₁::Mesh, mesh₂::Mesh) = merge(mesh₁, mesh₂)

function vizmany!(plot, objs)
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]
  pointsize = plot[:pointsize]
  segmentsize = plot[:segmentsize]

  object = Makie.@lift reduce(concat, $objs)
  colors = Makie.@lift mayberepeat($color, $objs)
  alphas = Makie.@lift mayberepeat($alpha, $objs)

  viz!(plot, object, color=colors, alpha=alphas, colorscheme=colorscheme, pointsize=pointsize, segmentsize=segmentsize)
end
