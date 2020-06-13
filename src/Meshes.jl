module Meshes

using StaticArrays: MVector # default buffer for coordinates

include("core.jl")
include("basic.jl")
include("geotraits.jl")

COMPILE_TIME_TRAITS = [:ismesh, :ndims,
                       :coordtype, :coordbuff,
                       :isstructured, :iscartesian,
                       :isrectilinear, :isregular]

# default versions for mesh instances
for TRAIT in COMPILE_TIME_TRAITS
  @eval $TRAIT(::M) where M = $TRAIT(M)
end

end # module
