# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

asalphas(alphas, _) = alphas
asalphas(::Nothing, values) = Colorfy.defaultalphas(values)

ascolorscheme(colorscheme, _) = colorscheme
ascolorscheme(::Nothing, values) = Colorfy.defaultscheme(values)

ascolorrange(colorrange, _) = colorrange
ascolorrange(::Nothing, _) = :extrema

# --------------------------------
# PROCESS COLORS PROVIDED BY USER
# --------------------------------

# convert user input to colors
function process(values::AbstractVector, colorscheme, colorrange, alphas)
  colorfy(
    values,
    alphas=asalphas(alphas, values),
    colorscheme=ascolorscheme(colorscheme, values),
    colorrange=ascolorrange(colorrange, values)
  )
end

process(value, colorscheme, colorrange, alphas) = process([value], colorscheme, colorrange, alphas) |> first
