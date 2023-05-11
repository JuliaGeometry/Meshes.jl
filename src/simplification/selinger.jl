# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Selinger()

Simplify geometries with Selinger algorithm, which attemps to
minimize the number of vertices and the deviation of vertices
to the resulting segments.

## References

* Selinger, P. 2003. [Potrace: A polygon-based tracing algorithm]
  (https://potrace.sourceforge.net/potrace.pdf)
"""
struct Selinger{T} <: SimplificationMethod
  ϵ::T
end

function simplify(chain::Chain{Dim,T}, method::Selinger) where {Dim,T}
  # retrieve parameters
  ϵ = method.ϵ

  # vertices as circular vector
  v = vertices(chain)
  p = isclosed(chain) ? v : CircularVector(v)

  # penalty for each possible segment
  n = length(p)
  P = Dict{Tuple{Int,Int},T}()
  for i in 1:n, o in 1:(n - 2)
    j = i + o
    i₊ = i + 1
    j₋ = j - 1
    jₙ = mod1(j, n)
    l = Line(p[i], p[j])
    δ = [evaluate(Euclidean(), p[k], l) for k in i₊:j₋]
    if all(<(ϵ), δ)
      dᵢⱼ = norm(p[j] - p[i])
      σᵢⱼ = o == 1 ? zero(T) : sqrt(sum(δ .^ 2) / length(δ))
      P[(i, jₙ)] = dᵢⱼ * σᵢⱼ
    end
  end

  # incidence matrix of directed graph
  I = spzeros(Int, n, n)
  for ij in keys(P)
    I[ij...] = 1
  end

  # shortest path with minimum penalty
  bestpath = dijkstra(I, 2, 1)
  bestlen = length(bestpath)
  bestpen = penalty(P, bestpath)
  for i in 2:(n - 1)
    path = dijkstra(I, i + 1, i)
    len = length(path)
    if len ≤ bestlen
      pen = penalty(P, path)
      if pen < bestpen
        bestpath = path
        bestlen = len
      end
    end
  end

  Chain(collect(v[bestpath]))
end

function dijkstra(I, s, t)
  score = Dict(s => 0)
  frontier = Dict(s => 0)
  parentof = Dict(s => 0)

  while !isempty(frontier)
    f, i = findmin(frontier)
    pop!(frontier, i)
    i == t && break
    for j in findall(I[i, :] .== 1)
      if j ∉ keys(score) || f + 1 < score[j]
        score[j] = f + 1
        frontier[j] = f + 1
        parentof[j] = i
      end
    end
  end

  buildcycle(parentof, t)
end

function buildcycle(parentof, t)
  p = buildpath(parentof, t)
  [p; first(p)]
end

function buildpath(parentof, j)
  p = [j]
  i = parentof[j]
  while i ≠ 0
    push!(p, i)
    i = parentof[i]
  end
  reverse(p)
end

penalty(P, path) = sum(P[(path[k], path[k + 1])] for k in 1:(length(path) - 1))
