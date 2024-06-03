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

function refine(grid::CartesianGrid{Dim}, method::RegularRefinement) where {Dim}
  factors = fitdims(method.factors, Dim)
  CartesianGrid(minimum(grid), maximum(grid), dims=size(grid) .* factors)
end

function refine(grid::RectilinearGrid{Datum,Dim}, method::RegularRefinement) where {Datum,Dim}
  factors = fitdims(method.factors, Dim)
  xyzₛ = xyz(grid)
  xyzₜ = ntuple(i -> _refinedims(xyzₛ[i], factors[i]), Dim)
  RectilinearGrid{Datum}(xyzₜ)
end

function refine(grid::StructuredGrid{Datum,Dim}, method::RegularRefinement) where {Datum,Dim}
  factors = fitdims(method.factors, Dim)
  XYZ′ = _XYZ(grid, factors)
  StructuredGrid{Datum}(XYZ′)
end

function _refinedims(x, f)
  x′ = mapreduce(vcat, 1:(length(x) - 1)) do i
    range(x[i], x[i + 1], f + 1)[begin:(end - 1)]
  end
  push!(x′, last(x))
  x′
end

function _XYZ(grid::StructuredGrid{Datum,2}, factors::Dims{2}) where {Datum}
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

function _XYZ(grid::StructuredGrid{Datum,3}, factors::Dims{3}) where {Datum}
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
