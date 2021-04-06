# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Data

A domain together with a table of values for each element of the domain.
The i-th row of the table is a vector of features associated with the
i-th element of the domain. If the domain is a mesh, then the table
stores all the properties of the elements.
"""
abstract type Data end

"""
    domain(data)

Return underlying domain of the `data`.
"""
function domain end

"""
    values(data)

Return the values of `data` as a table.
"""
values(data::Data)

"""
    constructor(D)

Return the constructor of the data type `D` as a function.
The function takes a `domain` and `table` as inputs and
combines them into an instance of the data type.
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
  cols = map(vars) do var
    var => Tables.getcolumn(ctable, var)
  end
  table = (; cols...)

  # combine the two with constructor
  constructor(D)(domain, table)
end

==(data₁::Data, data₂::Data) =
  domain(data₁) == domain(data₂) && values(data₁) == values(data₂)

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

# -------------------
# VARIABLE INTERFACE
# -------------------

"""
    getindex(data, var)

Returns the data for the variable `var` in `data` as a column vector.
"""
Base.getindex(data::Data, var::Symbol) =
  Tables.getcolumn(values(data), var)

Base.getindex(data::Data, var::String) =
  getindex(data, Symbol(var))

function variables(data::Data)
  s = Tables.schema(values(data))
  @. Variable(s.names, nonmissingtype(s.types))
end

# ----------
# UTILITIES
# ----------

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
  nelm = nelements(domain(data))
  Dim  = embeddim(domain(data))
  T    = coordtype(domain(data))
  print(io, "$nelm $name{$Dim,$T}")
end

function Base.show(io::IO, ::MIME"text/plain", data::Data)
  𝒟 = domain(data)
  𝒯 = values(data)
  s = Tables.schema(𝒯)
  vars = zip(s.names, s.types)
  println(io, data)
  println(io, "  variables")
  varlines = ["    └─$var ($V)" for (var,V) in vars]
  println(io, join(sort(varlines), "\n"))
  print(  io, "  domain: ", 𝒟)
end
