# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# auxiliary types for dispatch purposes
const GeometryOrDomain = Union{Geometry,Domain}
const CartesianOrProjected = Union{Cartesian,Projected}

include("utils/basic.jl")
include("utils/assert.jl")
include("utils/cmp.jl")
include("utils/units.jl")
include("utils/crs.jl")
include("utils/misc.jl")
