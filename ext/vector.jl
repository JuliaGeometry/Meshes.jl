# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{AbstractVector{Vec{Dim,ℒ}}}}) where {Dim,ℒ}
  vecs = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]
  segmentsize = plot[:segmentsize]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colormap, $colorrange, $alpha)

  # visualize as built-in arrows
  T = Unitful.numtype(ℒ)
  orig = Makie.@lift fill(zero(Makie.Point{Dim,T}), length($vecs))
  dirs = Makie.@lift asmakie.($vecs)
  if Dim == 2
    tipwidth = Makie.@lift 5 * $segmentsize
    shaftwidth = Makie.@lift 0.2 * $tipwidth
    Makie.arrows2d!(plot, orig, dirs, color=colorant, tipwidth=tipwidth, shaftwidth=shaftwidth)
  elseif Dim == 3
    tipradius = Makie.@lift 0.05 * $segmentsize
    shaftradius = Makie.@lift 0.5 * $tipradius
    Makie.arrows3d!(plot, orig, dirs, color=colorant, tipradius=tipradius, shaftradius=shaftradius)
  else
    error("not implemented")
  end
end
