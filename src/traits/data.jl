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
function constructor(::Type{Data}) end

# ----------
# FALLBACKS
# ----------

# implement Domain traits for convenience
embeddim(data::Data) = embeddim(domain(data))
coordtype(data::Data) = coordtype(domain(data))
nelements(data::Data) = nelements(domain(data))
coordinates!(buff, data::Data, ind::Int) =
  coordinates!(buff, domain(data), ind)

==(data₁::Data, data₂::Data) =
  domain(data₁) == domain(data₂) && values(data₁) == values(data₂)

# -----------------
# TABLES INTERFACE
# -----------------

Tables.istable(::Type{<:Data}) = true

Tables.rowaccess(data::Data) = true

function Tables.rows(data::Data)
  rows = Tables.rows(values(data))
  elms = domain(data)
  ((row..., elms[i]) for (i, row) in Iterators.enumerate(rows))
end

function Tables.schema(data::Data)
  geomtype = eltype(domain(data))
  schema = Tables.schema(values(data))
  names, types = schema.names, schema.types
  Tables.Schema((names..., :geometry), (types..., geomtype))
end

function Tables.materializer(::D) where {D<:Data}
  function materializer(stable)
    # build domain from geometry column
    elms = Tables.getcolumn(stable, :geometry)
    domain = GeometrySet(elms)

    # build table from remaining columns
    vars = setdiff(Tables.columnnames(stable), :geometry)
    cols = map(vars) do var
      var => Tables.getcolumn(stable, var)
    end
    table = (; cols...)

    # combine the two with constructor
    constructor(D)(domain, table)
  end
end

# -------------------
# VARIABLE INTERFACE
# -------------------

"""
    getindex(data, var)

Returns the data for the variable `var` in `data` as a column vector.
"""
Base.getindex(data::Data, var::Symbol) =
  Tables.getcolumn(values(data), var)

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
