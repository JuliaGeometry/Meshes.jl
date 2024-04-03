# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# preprocess colors provided by user
function process(values::AbstractVector, colorscheme, colorrange, alphas)
  valphas = isnothing(alphas) ? Colorfy.defaultalphas(values) : alphas
  vcolorscheme = isnothing(colorscheme) ? Colorfy.defaultscheme(values) : colorscheme
  vcolorrange = isnothing(colorrange) ? :extrema : colorrange
  colorfy(values, alphas=valphas, colorscheme=vcolorscheme, colorrange=vcolorrange)
end

process(value, colorscheme, colorrange, alphas) = process([value], colorscheme, colorrange, alphas) |> first
