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

# ----------
# FALLBACKS
# ----------

"""
    data₁ == data₂

Tells whether or not the `data₁` and `data₂` are equal.
"""
==(data₁::Data, data₂::Data) =
  domain(data₁) == domain(data₂) && values(data₁) == values(data₂)

# -----------
# TABLES API
# -----------

Tables.istable(::Type{<:Data}) = true
Tables.rowaccess(data::Data) = true
function Tables.schema(data::Data)
  geomtype = eltype(domain(data))
  schema = Tables.schema(values(data))
  names, types = schema.names, schema.types
  Tables.Schema((names..., :geometry), (types..., geomtype))
end
# Tables.columns(data::Data) = Tables.columns(values(data))
# Tables.columnnames(data::Data) = Tables.columnnames(values(data))
# Tables.getcolumn(data::Data, c::Symbol) = Tables.getcolumn(values(data), c)
# Tables.rows(data::Data) = Tables.rows(values(data))

# -------------
# VARIABLE API
# -------------

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

MIMES = [MIME"text/plain", MIME"text/html"]

for MIME in MIMES
  @eval function Base.show(io::IO, ::$MIME, data::Data)
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
end
