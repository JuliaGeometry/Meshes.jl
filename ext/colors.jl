# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# preprocess colors provided by user
process(values::AbstractVector, colorscheme, colorrange, alphas) = colorfy(values; alphas, colorscheme, colorrange)

process(value, colorscheme, colorrange, alphas) = process([value], colorscheme, colorrange, alphas) |> first
