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

Collect iterable `iter` at indices `inds` efficiently.
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

collectat(vec::AbstractVector, inds) = vec[inds]

"""
    XYZ(xyz)

Generate the coordinate arrays `XYZ` from the coordinate vectors `xyz`.
"""
@generated function XYZ(xyz::NTuple{Dim,AbstractVector}) where {Dim}
  exprs = ntuple(Dim) do d
    quote
      a = xyz[$d]
      A = Array{eltype(a),Dim}(undef, length.(xyz))
      @nloops $Dim i A begin
        @nref($Dim, A, i) = a[$(Symbol(:i_, d))]
      end
      A
    end
  end
  Expr(:tuple, exprs...)
end

"""
    roundcoords(point, r=RoundNearest; digits=0, base=10)
    roundcoords(point, r=RoundNearest; sigdigits=0)

Rounds the coordinates of a `point` to specified presicion.
"""
function roundcoords(p::Point, r::RoundingMode=RoundNearest; kwargs...)
  c = coords(p)
  x = CoordRefSystems.values(c)
  x′ = round.(eltype(x), x, r; kwargs...)
  c′ = CoordRefSystems.constructor(c)(x′...)
  Point(c′)
end
