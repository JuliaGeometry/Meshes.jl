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

  Dim ∈ (2, 3) || error("not implemented")

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colormap, $colorrange, $alpha)

  # visualize as built-in arrows
  T = Unitful.numtype(ℒ)
  orig = Makie.@lift fill(zero(Makie.Point{Dim,T}), length($vecs))
  dirs = Makie.@lift asmakie.($vecs)
  size = Makie.@lift 0.1 * $segmentsize
  Makie.arrows!(plot, orig, dirs, color=colorant, arrowsize=size)
end
