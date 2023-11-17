# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module MeshesMakieExt

using Meshes

using Makie: cgrad
using Makie: coloralpha
using Makie.Colors: Colorant

import Meshes: viz, viz!
import Meshes: ascolors
import Meshes: defaultscheme
import Makie

Makie.@recipe(Viz, object) do scene
  Makie.Attributes(
    color=:slategray3,
    alpha=nothing,
    colorscheme=nothing,
    facetcolor=:gray30,
    showfacets=false,
    pointsize=2,
    segmentsize=1.5
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
include("fallbacks.jl")
include("optimized.jl")

end
