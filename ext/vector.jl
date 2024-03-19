# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{AbstractVector{Vec{Dim,T}}}}) where {Dim,T}
  vecs = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]
  segmentsize = plot[:segmentsize]

  Dim ∈ (2, 3) || error("not implemented")

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colorscheme, $alpha)

  # visualize as built-in arrows
  orig = Makie.@lift fill(zero(Makie.Point{Dim,T}), length($vecs))
  dirs = Makie.@lift asmakie.($vecs)
  size = Makie.@lift 0.1 * $segmentsize
  Makie.arrows!(plot, orig, dirs, color=colorant, arrowsize=size)
end
