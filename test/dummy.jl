# dummy type implementing the Data trait
struct DummyData{𝒟,𝒱} <: Data
  domain::𝒟
  values::𝒱
end

Meshes.domain(data::DummyData) = getfield(data, :domain)

function Meshes.values(data::DummyData, rank=nothing)
  domain = getfield(data, :domain)
  values = getfield(data, :values)
  r = isnothing(rank) ? paramdim(domain) : rank
  haskey(values, r) ? values[r] : nothing
end

Meshes.constructor(::Type{D}) where {D<:DummyData} = DummyData
