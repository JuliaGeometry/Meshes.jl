# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    constructor(G)

Given a (parametric) type `G{T₁,T₂,...}`, return the type `G`.
"""
constructor(G::Type) = getfield(Meshes, nameof(G))

"""
    fitdims(dims, D)

Fit tuple `dims` to a given length `D` by repeating the last dimension.
"""
function fitdims(dims::Dims{N}, D) where {N}
  ntuple(i -> i ≤ N ? dims[i] : last(dims), D)
end

"""
    collectat(iter, inds)

Collect iterator `iter` at indices `inds` without materialization.
"""
function collectat(iter, inds)
  if isempty(inds)
    eltype(iter)[]
  else
    m = maximum(inds)
    e = Iterators.enumerate(iter)
    w = Iterators.takewhile(x -> (first(x) ≤ m), e)
    f = Iterators.filter(x -> (first(x) ∈ inds), w)
    map(last, f)
  end
end

"""
    XYZ(xyz)

Generate the coordinate arrays `XYZ` from the coordinate vectors `xyz`.
"""
@generated function XYZ(xyz::NTuple{Dim,<:AbstractVector{T}}) where {Dim,T}
  exprs = ntuple(Dim) do d
    quote
      a = xyz[$d]
      A = Array{T,Dim}(undef, length.(xyz))
      @nloops $Dim i A begin
        @nref($Dim, A, i) = a[$(Symbol(:i_, d))]
      end
      A
    end
  end
  Expr(:tuple, exprs...)
end
