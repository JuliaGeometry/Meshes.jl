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
    vals₁ = values(data₁, rank)
    vals₂ = values(data₂, rank)
    if !isequal(vals₁, vals₂)
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
Tables.columnaccess(::Type{<:Data}) = true
Tables.columns(data::Data) = DataCols(data)
Tables.rowaccess(::Type{<:Data}) = true
Tables.rows(data::Data) = DataRows(DataCols(data))
Tables.schema(data::Data) = Tables.schema(DataCols(data))
Tables.materializer(::Type{D}) where {D<:Data} = constructor(D)

# wrapper type for cols of the data table
struct DataCols{D<:Data,C,G}
  tcols::C
  ncols::Int
  names::Vector{Symbol}
  domain::Vector{G}

  function DataCols(data::D) where {D<:Data}
    tcols = Tables.columns(values(data))
    names = [Tables.columnnames(tcols)..., :geometry]
    dom   = collect(domain(data))
    new{D,typeof(tcols),eltype(dom)}(tcols, length(names), names, dom)
  end
end

Tables.istable(::Type{<:DataCols}) = true
Tables.columnaccess(::Type{<:DataCols}) = true
Tables.columns(cols::DataCols) = cols
Tables.rowaccess(::Type{<:DataCols}) = true
Tables.rows(cols::DataCols) = DataRows(cols)
Tables.columnnames(cols::DataCols) = cols.names

function Tables.getcolumn(cols::DataCols, i::Int)
  1 ≤ i ≤ cols.ncols || error("Table has no column with index $i.")
  Tables.getcolumn(cols, cols.names[i])
end

function Tables.getcolumn(cols::DataCols, nm::Symbol)
  nm ∉ cols.names && error("Table has no column $nm.")
  nm == :geometry ? cols.domain : Tables.getcolumn(cols.tcols, nm)
end

function Tables.schema(cols::DataCols)
  types   = Tables.schema(cols.tcols).types
  geotype = eltype(cols.domain)
  Tables.Schema(cols.names, [types..., geotype])
end

Tables.materializer(::Type{DataCols{D,C,G}}) where {D<:Data,C,G} =
  constructor(D)

# wrapper type for row of the data table
struct DataRow{C<:DataCols}
  cols::C
  ind::Int
  ncols::Int

  function DataRow(cols::C, ind::Int) where {C<:DataCols}
    new{C}(cols, ind, cols.ncols)
  end
end

# Iteration interface
Base.iterate(row::DataRow, state::Int=1) =
  state > row.ncols ? nothing : (row[state], state + 1)

Base.length(row::DataRow) = row.ncols
Base.IteratorSize(::Type{<:DataRow}) = Base.HasLength()
Base.IteratorEltype(::Type{<:DataRow}) = Base.EltypeUnknown()

# Indexing interface
Base.firstindex(::DataRow) = 1
Base.lastindex(row::DataRow) = row.ncols
Base.eachindex(row::DataRow) = 1:row.ncols
Base.getindex(row::DataRow, i::Int) = Tables.getcolumn(row, i)

# Tables.jl row interface
Tables.columnnames(row::DataRow) = Tables.columnnames(row.cols)
Tables.getcolumn(row::DataRow, i::Int) =
  Tables.getcolumn(row.cols, i)[row.ind]
Tables.getcolumn(row::DataRow, nm::Symbol) =
  Tables.getcolumn(row.cols, nm)[row.ind]

# wrapper type for rows of the data table
struct DataRows{C<:DataCols}
  cols::C
  nrows::Int

  function DataRows(cols::C) where {C<:DataCols}
    new{C}(cols, length(cols.domain))
  end
end

# Iteration interface
Base.iterate(rows::DataRows, state::Int=1) =
  state > rows.nrows ? nothing : (rows[state], state + 1)

Base.length(rows::DataRows) = rows.nrows
Base.eltype(::Type{DataRows{T}}) where {T} = DataRow{T}
Base.IteratorSize(::Type{<:DataRows}) = Base.HasLength()
Base.IteratorEltype(::Type{<:DataRows}) = Base.HasEltype()

# Indexing interface
Base.firstindex(::DataRows) = 1
Base.lastindex(rows::DataRows) = rows.nrows
Base.eachindex(rows::DataRows) = 1:rows.nrows
Base.getindex(rows::DataRows, i::Int) = DataRow(rows.cols, i)

# Tables.jl interface
Tables.isrowtable(::Type{<:DataRows}) = true
Tables.columnaccess(::Type{<:DataRows}) = true
Tables.columns(rows::DataRows) = Tables.columns(rows.cols)
Tables.columnnames(rows::DataRows) = Tables.columnnames(rows.cols)
Tables.getcolumn(rows::DataRows, i::Int) = Tables.getcolumn(rows.cols, i)
Tables.getcolumn(rows::DataRows, nm::Symbol) = Tables.getcolumn(rows.cols, nm)
Tables.materializer(::Type{DataRows{C}}) where {C<:DataCols} = Tables.materializer(C)
Tables.schema(rows::DataRows) = Tables.schema(rows.cols)

# data table is compatible with the Queryverse
TableTraits.isiterabletable(data::Data) = true
IIE.getiterator(data::Data) = Tables.datavaluerows(Tables.rows(data))
IIE.isiterable(data::Data) = true

# -----------------
# COLUMN INTERFACE
# -----------------

function Base.getindex(data::Data, col::Symbol)
  if col == :geometry
    collect(domain(data))
  else
    cols = Tables.columns(values(data))
    Tables.getcolumn(cols, col)
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
      sche = Tables.schema(𝒯)
      vars = zip(sche.names, sche.types)
      println(io, "  variables (rank $rank)")
      varlines = ["    └─$var ($V)" for (var,V) in vars]
      println(io, join(sort(varlines), "\n"))
    end
  end
  print(io, "  domain: ", 𝒟)
end
