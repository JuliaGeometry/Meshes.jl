# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    View(discretization, indices)

A partial view of a `discretization` containing only the elements at `indices`.
"""
struct View{Dim,T,D<:Discretization{Dim,T},I} <: Discretization{Dim,T}
  disc::D
  inds::I
end

# convenience functions
Base.view(disc::Discretization, inds) = View(disc, inds)

# -------------------------
# DISCRETIZATION INTERFACE
# -------------------------

Base.getindex(v::View, ind::Int) = getindex(v.disc, v.inds[ind])

nelements(v::View) = length(v.inds)

coordinates!(buff, v::View, ind::Int) = 
  coordinates!(buff, v.disc, v.inds[ind])

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, v::View)
  disc  = v.disc
  nelms = length(v.inds)
  print(io, "$nelms View{$disc}")
end
