# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

include("traits/variable.jl")
include("traits/domain.jl")
include("traits/data.jl")
include("traits/optional.jl")

# type alias for convenience
const DomainOrData = Union{Domain,Data}
