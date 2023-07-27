# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module MeshesMakieExt

using Meshes

using Makie.Colors: Colorant
using Makie.Colors: protanopic, coloralpha
using Makie.Colors: distinguishable_colors
using Makie.ColorSchemes: colorschemes

import Meshes: viz, viz!
import Makie

@Makie.recipe(Viz, object) do scene
  Makie.Attributes(
    color       = :slategray3,
    alpha       = 1.0,
    colorscheme = nothing,
    facetcolor  = :gray30,
    showfacets  = false,
    pointsize   = 12,
    segmentsize = 1.5
  )
end

# color handling
include("colors.jl")

# utilities
include("utils.jl")

# viz recipes
include("simplemesh.jl")
include("cartesiangrid.jl")
include("geometryset.jl")
include("partition.jl")
include("fallbacks.jl")
include("optimized.jl")

end