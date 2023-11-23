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
    pointsize=2,
    segmentsize=1.5,
    showfacets=false,
    facetcolor=:gray30
  )
end

# choose between 2D and 3D axis
Makie.args_preferred_axis(::Geometry{Dim}) where {Dim} = Dim === 2 ? Makie.Axis : Makie.LScene
Makie.args_preferred_axis(::Domain{Dim}) where {Dim} = Dim === 2 ? Makie.Axis : Makie.LScene

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
