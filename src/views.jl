# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# -------------
# DOMAIN VIEWS
# -------------

"""
    DomainView(domain, indices)

A partial view of a `domain` containing only the elements at `indices`.
"""
struct DomainView{Dim,T,D<:Domain{Dim,T},I} <: Domain{Dim,T}
  domain::D
  inds::I
end

# -----------------
# DOMAIN INTERFACE
# -----------------

Base.getindex(v::DomainView, ind::Int) =
  getindex(v.domain, v.inds[ind])

nelements(v::DomainView) = length(v.inds)

coordinates!(buff, v::DomainView, ind::Int) = 
  coordinates!(buff, v.domain, v.inds[ind])

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, v::DomainView)
  domain  = v.domain
  nelms = length(v.inds)
  print(io, "$nelms View{$domain}")
end

# -----------
# DATA VIEWS
# -----------

# helper function for table view
function viewtable(table, rows, cols)
  t = Tables.columns(table)
  v = map(cols) do c
    col = Tables.getcolumn(t, c)
    c => view(col, rows)
  end
  (; v...)
end

"""
    DataView(data, inds, vars)

Return a view of `data` at indices `inds` and variables `vars`.
"""
struct DataView{D<:Data,I,V} <: Data
  data::D
  inds::I
  vars::V
end

# ---------------
# DATA INTERFACE
# ---------------

domain(v::DataView) = view(domain(v.data), v.inds)
values(v::DataView) = viewtable(values(v.data), v.inds, v.vars)
function constructor(::Type{DataView{D,I,V}}) where {D<:Data,I,V}
  function ctor(domain, table)
    data = constructor(D)(domain, table)
    inds = 1:nelements(domain)
    vars = Tables.schema(table).names
    DataView(data, inds, vars)
  end
end

# specialize methods for performance
==(v₁::DataView, v₂::DataView) =
  v₁.data == v₂.data && v₁.inds == v₂.inds && v₁.vars == v₂.vars

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, v::DataView)
  data  = v.data
  nelms = length(v.inds)
  print(io, "$nelms View{$data}")
end
