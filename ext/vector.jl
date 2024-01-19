# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{AbstractVector{Vec{Dim,T}}}}) where {Dim,T}
  vecs = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]

  if Dim ∉ (2, 3)
    error("not implemented")
  end

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colorscheme, $alpha)

  # visualize as built-in arrows
  origins = Makie.@lift fill(zero(Makie.Point{Dim,T}), length($vecs))
  directions = Makie.@lift asmakie.($vecs)
  Makie.arrows!(plot, origins, directions, color=colorant)
end
