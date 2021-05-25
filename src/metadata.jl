# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Metadata(domain, values)

A `domain` together with a dictionary of data `values`. For each rank `r`
(or parametric dimension) there can exist a corresponding Tables.jl table
`values[r]`. The helper function [`metadata`](@ref) is recommended instead
of the raw constructor of the type.
"""
struct Metadata{D<:Domain,V<:Dict} <: Data
  domain::D
  values::V
end

domain(data::Metadata) = getfield(data, :domain)

function values(data::Metadata, rank=nothing)
  domain = getfield(data, :domain)
  values = getfield(data, :values)
  r = isnothing(rank) ? paramdim(domain) : rank
  haskey(values, r) ? values[r] : nothing
end

constructor(::Type{D}) where {D<:Metadata} = Metadata

# ----------------
# HELPER FUNCTION
# ----------------

"""
    metadata(domain, values)

Create spatial data from a `domain` implementing the
the [`Domain`](@ref) trait and a dictionary of data
`values` where `values[r]` holds a Tables.jl table
for the rank `r`.
"""
metadata(domain::Domain, values::Dict) = Metadata(domain, values)

"""
    metadata(vertices, elements, values)

Create spatial data from a [`SimpleMesh`](@ref) with
`vertices` and `elements`, and a dictionary of data
`values`.
"""
metadata(vertices, elements, values) =
  metadata(SimpleMesh(vertices, elements), values)

"""
    metadata(domain; vtable, etable)

Create spatial data from a `domain`, a table `vtable`
with data for the vertices, and a table `etable` with
data for the elements.

## Example

```julia
metadata(CartesianGrid(10, 10),
  etable = (temperature=rand(100), pressure=rand(100))
)
```
"""
function metadata(domain; vtable=nothing, etable=nothing)
  d = paramdim(domain)
  values = if !isnothing(vtable) && !isnothing(etable)
    Dict(0 => vtable, d => etable)
  elseif isnothing(vtable)
    Dict(d => etable)
  elseif isnothing(etable)
    Dict(0 => vtable)
  else
    throw(ArgumentError("missing data tables"))
  end
  metadata(domain, values)
end

"""
    metadata(vertices, elements; vtable, etable)

Create spatial data from a [`SimpleMesh`](@ref) with
`vertices` and `elements`, a table `vtable` with data
for the vertices, and a table `etable` with data for
the elements.

## Example

```julia
# vertices and elements
vertices = Point2[(0,0),(1,0),(1,1),(0,1)]
elements = connect.([(1,2,3),(3,4,1)])

# attach data to mesh
metadata(vertices, elements,
  vtable = (temperature=[1.0,2.0,3.0,4.0], pressure=[4.0,3.0,2.0,1.0]),
  etable = (quality=["A","B"], state=[true,false])
)
```
"""
metadata(vertices, elements; vtable=nothing, etable=nothing) =
  metadata(SimpleMesh(vertices, elements); vtable=vtable, etable=etable)
