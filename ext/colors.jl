# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# type alias to reduce typing
const V{T} = AbstractVector{<:T}

# convert value to colorant, optionally using color scheme object
ascolors(values::V{Symbol}, scheme) = ascolors(string.(values), scheme)
ascolors(values::V{AbstractString}, scheme) = parse.(Ref(Colorant), values)
ascolors(values::V{Number}, scheme) = get(scheme, values, :extrema)
ascolors(values::V{Colorant}, scheme) = values

# convert color scheme name to color scheme object
ascolorscheme(name::Symbol) = colorschemes[name]
ascolorscheme(name::AbstractString) = ascolorscheme(Symbol(name))
ascolorscheme(scheme) = scheme

# default colorscheme for vector of values
defaultscheme(values) = colorschemes[:viridis]

# --------------------------------
# PROCESS COLORS PROVIDED BY USER
# --------------------------------

# STEP 1: convert user input to colors
function tocolors(values, scheme)
  # find invalid and valid indices
  iinds = findall(ismissing, values)
  vinds = setdiff(1:length(values), iinds)

  # invalid values are assigned full transparency
  icolors = parse(Colorant, "rgba(0,0,0,0)")

  # valid values are assigned colors from scheme
  vals = coalesce.(values[vinds])
  vscheme = isnothing(scheme) ? defaultscheme(vals) : ascolorscheme(scheme)
  vcolors = ascolors(vals, vscheme)

  # build final vector of colors
  colors = Vector{Colorant}(undef, length(values))
  colors[iinds] .= icolors
  colors[vinds] .= vcolors

  colors
end

# STEP 2: add transparency to colors
setalpha(colors, alphas) = coloralpha.(colors, alphas)
setalpha(colors, alpha::Number) = coloralpha.(colors, alpha)
setalpha(colors, ::Nothing) = colors

process(value, scheme, alphas) = process([value], scheme, alphas) |> first
process(values::V, scheme, alphas) = setalpha(tocolors(values, scheme), alphas)
