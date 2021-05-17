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

element(v::DomainView, ind::Int) =
  element(v.domain, v.inds[ind])

nelements(v::DomainView) = length(v.inds)

centroid(v::DomainView, ind::Int) = 
  centroid(v.domain, v.inds[ind])

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

domain(v::DataView) = view(domain(getfield(v, :data)),
                           getfield(v, :inds))

values(v::DataView) = viewtable(values(getfield(v, :data)),
                                getfield(v, :inds),
                                getfield(v, :vars))

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
  getfield(v₁, :data) == getfield(v₂, :data) &&
  getfield(v₁, :inds) == getfield(v₂, :inds) &&
  getfield(v₁, :vars) == getfield(v₂, :vars)

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, v::DataView)
  data  = getfield(v, :data)
  nelms = length(getfield(v, :inds))
  print(io, "$nelms View{$data}")
end
