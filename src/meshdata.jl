# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MeshData(domain, values)

A `domain` together with a dictionary of data `values`. For each rank `r`
(or parametric dimension) there can exist a corresponding Tables.jl table
`values[r]`. The helper function [`meshdata`](@ref) is recommended instead
of the raw constructor of the type.
"""
struct MeshData{D<:Domain,V<:Dict} <: Data
  domain::D
  values::V
end

domain(data::MeshData) = getfield(data, :domain)

function values(data::MeshData, rank=nothing)
  domain = getfield(data, :domain)
  values = getfield(data, :values)
  r = isnothing(rank) ? paramdim(domain) : rank
  haskey(values, r) ? values[r] : nothing
end

constructor(::Type{D}) where {D<:MeshData} = MeshData

# ----------------
# HELPER FUNCTION
# ----------------

"""
    meshdata(domain, values)

Create spatial data from a `domain` implementing the
the [`Domain`](@ref) trait and a dictionary of data
`values` where `values[r]` holds a Tables.jl table
for the rank `r`.

## Example

```julia
# attach temperature and pressure to grid elements
meshdata(CartesianGrid(10,10),
  Dict(2 => (temperature=rand(100), pressure=rand(100)))
)
```
"""
meshdata(domain::Domain, values::Dict) = MeshData(domain, values)

"""
    meshdata(vertices, elements, values)

Create spatial data from a [`SimpleMesh`](@ref) with
`vertices` and `elements`, and a dictionary of data
`values`.

## Example

```julia
# vertices and elements
vertices = Point2[(0,0),(1,0),(1,1),(0,1)]
elements = connect.([(1,2,3),(3,4,1)])

# attach data to ranks 0 and 2
meshdata(vertices, elements,
  Dict(
    0 => (temperature=[1.0,2.0,3.0,4.0], pressure=[4.0,3.0,2.0,1.0]),
    2 => (quality=["A","B"], state=[true,false])
  )
)
```
"""
meshdata(vertices, elements, values) =
  meshdata(SimpleMesh(vertices, elements), values)

"""
    meshdata(domain; vtable, etable)

Create spatial data from a `domain`, a table `vtable`
with data for the vertices, and a table `etable` with
data for the elements.

## Example

```julia
meshdata(CartesianGrid(10,10),
  etable = (temperature=rand(100), pressure=rand(100))
)
```
"""
function meshdata(domain; vtable=nothing, etable=nothing)
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
  meshdata(domain, values)
end

"""
    meshdata(vertices, elements; vtable, etable)

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
meshdata(vertices, elements,
  vtable = (temperature=[1.0,2.0,3.0,4.0], pressure=[4.0,3.0,2.0,1.0]),
  etable = (quality=["A","B"], state=[true,false])
)
```
"""
meshdata(vertices, elements; vtable=nothing, etable=nothing) =
  meshdata(SimpleMesh(vertices, elements); vtable=vtable, etable=etable)
