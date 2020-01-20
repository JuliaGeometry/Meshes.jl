module Meshes

using StaticArrays: MVector # default buffer for coordinates

include("core.jl")
include("basic.jl")
include("geoprops.jl")

COMPILE_TIME_TRAITS = [:ismesh, :ndims, :coordtype, :coordbuff,
                       :isstructured, :isregular]

# default versions for mesh instances
for TRAIT in COMPILE_TIME_TRAITS
  @eval $TRAIT(::M) where M = $TRAIT(M)
end

end # module
