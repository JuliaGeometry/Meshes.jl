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
ascolorscheme(name::Symbol) = cgrad(name)
ascolorscheme(name::AbstractString) = ascolorscheme(Symbol(name))
ascolorscheme(scheme) = scheme

# default color scheme for vector of values
defaultscheme(values) = cgrad(:viridis)

# add transparency to colors
setalpha(colors, alphas) = coloralpha.(colors, alphas)
setalpha(colors, ::Nothing) = colors

# --------------------------------
# PROCESS COLORS PROVIDED BY USER
# --------------------------------

# convert user input to colors
function process(values::V, scheme, alphas)
  # find invalid and valid indices
  isinvalid(v) = ismissing(v) || (v isa Number && isnan(v))
  iinds = findall(isinvalid, values)
  vinds = setdiff(1:length(values), iinds)

  # invalid values are assigned full transparency
  icolors = parse(Colorant, "rgba(0,0,0,0)")

  # valid values are assigned colors from scheme
  vals = coalesce.(values[vinds])
  vscheme = isnothing(scheme) ? defaultscheme(vals) : ascolorscheme(scheme)
  vcolors = setalpha(ascolors(vals, vscheme), alphas)

  # build final vector of colors
  colors = Vector{Colorant}(undef, length(values))
  colors[iinds] .= icolors
  colors[vinds] .= vcolors

  colors
end

process(value, scheme, alphas) = process([value], scheme, alphas) |> first
