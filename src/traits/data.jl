# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Data

A domain implementing the [`Domain`](@ref) trait together with tables
of values for geometries of the domain.

See also [`meshdata`](@ref).
"""
abstract type Data end

"""
    domain(data)

Return underlying domain of the `data`.
"""
function domain end

"""
    values(data, [rank])

Return the values of `data` for a given `rank` as a table.

The rank is a non-negative integer that specifies the
parametric dimension of the geometries of interest:

* 0 - points
* 1 - segments
* 2 - triangles, quadrangles, ...
* 3 - tetrahedrons, hexahedrons, ...

If the rank is not specified, it is assumed to be the rank
of the elements of the domain.
"""
function values end

"""
    constructor(D)

Return the constructor of the data type `D` as a function.
The function takes a domain and a dictionary of tables as
inputs and combines them into an instance of the data type.
"""
function constructor end

# ----------
# FALLBACKS
# ----------

function (D::Type{<:Data})(stable)
  # build domain from geometry column
  ctable = Tables.columns(stable)
  elms   = Tables.getcolumn(ctable, :geometry)
  domain = Collection(elms)

  # build table from remaining columns
  vars = setdiff(Tables.columnnames(ctable), (:geometry,))
  cols = [var => Tables.getcolumn(ctable, var) for var in vars]
  table = (; cols...)

  # data table for elements
  values = Dict(paramdim(domain) => table)

  # combine the two with constructor
  constructor(D)(domain, values)
end

function ==(data₁::Data, data₂::Data)
  # must have the same domain
  if domain(data₁) != domain(data₂)
    return false
  end

  # must have the same data tables
  for rank in 0:paramdim(domain(data₁))
    if values(data₁, rank) != values(data₂, rank)
      return false
    end
  end

  return true
end

# implement Domain traits for convenience
embeddim(data::Data) = embeddim(domain(data))
coordtype(data::Data) = coordtype(domain(data))
nelements(data::Data) = nelements(domain(data))
centroid(data::Data, ind) = centroid(domain(data), ind)

# -----------------
# TABLES INTERFACE
# -----------------

Tables.istable(::Type{<:Data}) = true

Tables.rowaccess(::Type{<:Data}) = true

Tables.rows(data::Data) = DataRows(domain(data), Tables.rows(values(data)))

# wrapper type for rows of the data table
# so that we can easily inform the schema
struct DataRows{𝒟,𝒯}
  domain::𝒟
  rtable::𝒯
end

function Base.iterate(rows::DataRows, state=1)
  if state > nelements(rows.domain)
    nothing
  else
    row, _ = iterate(rows.rtable, state)
    elm, _ = iterate(rows.domain, state)
    (; NamedTuple(row)..., geometry=elm), state + 1
  end
end

Base.length(rows::DataRows) = nelements(rows.domain)

function Tables.schema(rows::DataRows)
  geomtype = eltype(rows.domain)
  schema = Tables.schema(rows.rtable)
  names, types = schema.names, schema.types
  Tables.Schema((names..., :geometry), (types..., geomtype))
end

# data table is compatible with the Queryverse
TableTraits.isiterabletable(data::Data) = true
IIE.getiterator(data::Data) = Tables.datavaluerows(Tables.rows(data))
IIE.isiterable(data::Data) = true

Tables.materializer(D::Type{<:Data}) = D

# -----------------
# COLUMN INTERFACE
# -----------------

function Base.getindex(data::Data, col::Symbol)
  if col == :geometry
    collect(domain(data))
  else
    Tables.getcolumn(values(data), col)
  end
end

Base.getindex(data::Data, col::String) =
  getindex(data, Symbol(col))

Base.getproperty(data::Data, col::Symbol) =
  getindex(data, col)

# -------------------
# VARIABLE INTERFACE
# -------------------

"""
    variables(data)

Returns the variables stored in `data` as a vector of
[`Variable`](@ref).
"""
function variables(data::Data)
  s = Tables.schema(values(data))
  @. Variable(s.names, nonmissingtype(s.types))
end

"""
    asarray(data, var)

Returns the data for the variable `var` in `data` as a Julia array
with size equal to the size of the underlying domain if the size is
defined, otherwise returns a vector.
"""
function asarray(data::Data, var::Symbol)
  D = domain(data)
  hassize = hasmethod(size, (typeof(D),))
  hassize ? reshape(data[var], size(D)) : data[var]
end

asarray(data::Data, var::String) =
  asarray(data, Symbol(var))

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, data::Data)
  name = nameof(typeof(data))
  𝒟    = domain(data)
  n    = nelements(𝒟)
  Dim  = embeddim(𝒟)
  T    = coordtype(𝒟)
  print(io, "$n $name{$Dim,$T}")
end

function Base.show(io::IO, ::MIME"text/plain", data::Data)
  name = nameof(typeof(data))
  𝒟    = domain(data)
  n    = nelements(𝒟)
  Dim  = embeddim(𝒟)
  T    = coordtype(𝒟)
  println(io, "$n $name{$Dim,$T}")
  for rank in 0:paramdim(𝒟)
    𝒯 = values(data, rank)
    if !isnothing(𝒯)
      sche = Tables.schema(Tables.rows(𝒯))
      vars = zip(sche.names, sche.types)
      println(io, "  variables (rank $rank)")
      varlines = ["    └─$var ($V)" for (var,V) in vars]
      println(io, join(sort(varlines), "\n"))
    end
  end
  print(io, "  domain: ", 𝒟)
end
