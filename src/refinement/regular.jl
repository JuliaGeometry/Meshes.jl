# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RegularRefinement(f₁, f₂, ..., fₙ)

Refine each dimension of the grid by given factors `f₁`, `f₂`, ..., `fₙ`.

## Examples

```julia
refine(grid2D, RegularRefinement(2, 3))
refine(grid3D, RegularRefinement(2, 3, 1))
```
"""
struct RegularRefinement{N} <: RefinementMethod
  factors::Dims{N}
end

RegularRefinement(factors::Vararg{Int,N}) where {N} = RegularRefinement(factors)

function refine(grid::RegularGrid{<:𝔼,<:CartesianOrProjected}, method::RegularRefinement)
  factors = fitdims(method.factors, paramdim(grid))
  RegularGrid(minimum(grid), maximum(grid), dims=size(grid) .* factors)
end

function refine(grid::RectilinearGrid, method::RegularRefinement)
  factors = fitdims(method.factors, paramdim(grid))
  xyzₛ = xyz(grid)
  xyzₜ = ntuple(i -> _refinedims(xyzₛ[i], factors[i]), paramdim(grid))
  RectilinearGrid{manifold(grid),crs(grid)}(xyzₜ)
end

function refine(grid::StructuredGrid, method::RegularRefinement)
  factors = fitdims(method.factors, paramdim(grid))
  XYZ′ = _XYZ(grid, factors)
  StructuredGrid{manifold(grid),crs(grid)}(XYZ′)
end

refine(grid::TransformedGrid, method::RegularRefinement) =
  TransformedGrid(refine(parent(grid), method), transform(grid))

function _refinedims(x, f)
  x′ = mapreduce(vcat, 1:(length(x) - 1)) do i
    range(x[i], x[i + 1], f + 1)[begin:(end - 1)]
  end
  push!(x′, last(x))
  x′
end

_XYZ(grid::StructuredGrid, factors::Dims) = _XYZ(grid, Val(paramdim(grid)), factors)

function _XYZ(grid::StructuredGrid, ::Val{2}, factors::Dims{2})
  T = numtype(lentype(grid))
  fᵢ, fⱼ = factors
  sᵢ, sⱼ = size(grid)
  us = 0:T(1 / fᵢ):1
  vs = 0:T(1 / fⱼ):1
  catᵢ(A...) = cat(A..., dims=Val(1))
  catⱼ(A...) = cat(A..., dims=Val(2))

  mat(quad) = [to(quad(u, v)) for u in us, v in vs]
  M = [mat(grid[i, j]) for i in 1:sᵢ, j in 1:sⱼ]

  C = mapreduce(catⱼ, 1:sⱼ) do j
    Mⱼ = mapreduce(catᵢ, 1:sᵢ) do i
      Mᵢⱼ = M[i, j]
      i == sᵢ ? Mᵢⱼ : Mᵢⱼ[begin:(end - 1), :]
    end
    j == sⱼ ? Mⱼ : Mⱼ[:, begin:(end - 1)]
  end

  X = getindex.(C, 1)
  Y = getindex.(C, 2)

  (X, Y)
end

function _XYZ(grid::StructuredGrid, ::Val{3}, factors::Dims{3})
  T = numtype(lentype(grid))
  fᵢ, fⱼ, fₖ = factors
  sᵢ, sⱼ, sₖ = size(grid)
  us = 0:T(1 / fᵢ):1
  vs = 0:T(1 / fⱼ):1
  ws = 0:T(1 / fₖ):1
  catᵢ(A...) = cat(A..., dims=Val(1))
  catⱼ(A...) = cat(A..., dims=Val(2))
  catₖ(A...) = cat(A..., dims=Val(3))

  mat(hex) = [to(hex(u, v, w)) for u in us, v in vs, w in ws]
  M = [mat(grid[i, j, k]) for i in 1:sᵢ, j in 1:sⱼ, k in 1:sₖ]

  C = mapreduce(catₖ, 1:sₖ) do k
    Mₖ = mapreduce(catⱼ, 1:sⱼ) do j
      Mⱼₖ = mapreduce(catᵢ, 1:sᵢ) do i
        Mᵢⱼₖ = M[i, j, k]
        i == sᵢ ? Mᵢⱼₖ : Mᵢⱼₖ[begin:(end - 1), :, :]
      end
      j == sⱼ ? Mⱼₖ : Mⱼₖ[:, begin:(end - 1), :]
    end
    k == sₖ ? Mₖ : Mₖ[:, :, begin:(end - 1)]
  end

  X = getindex.(C, 1)
  Y = getindex.(C, 2)
  Z = getindex.(C, 3)

  (X, Y, Z)
end
