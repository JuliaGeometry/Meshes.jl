# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    laplacematrix(mesh; weights=:uniform)

The Laplace-Beltrami (a.k.a. Laplacian) matrix of the `mesh`.
Optionally specify the discretization `weights`.

## Weights

* `:uniform`   - `Lᵢⱼ = 1 / |𝒩(i)|, ∀j ∈ 𝒩(i)`
* `:cotangent` - `Lᵢⱼ = cot(αᵢⱼ) + cot(βᵢⱼ), ∀j ∈ 𝒩(i)`

## References

* Botsch et al. 2010. [Polygon Mesh Processing](http://www.pmp-book.org).

* Pinkall, U. & Polthier, K. 1993. [Computing discrete minimal surfaces and their conjugates]
  (https://projecteuclid.org/journals/experimental-mathematics/volume-2/issue-1/Computing-discrete-minimal-surfaces-and-their-conjugates/em/1062620735.full).
"""
function laplacematrix(mesh; weights=:uniform)
  # convert to half-edge topology
  ℳ = topoconvert(HalfEdgeTopology, mesh)

  # retrieve adjacency relation
  𝒩 = Adjacency{0}(topology(ℳ))

  # initialize matrix
  n = nvertices(ℳ)
  L = spzeros(n, n)

  # fill matrix with weights
  if weights == :uniform
    uniformlaplacian!(L, 𝒩)
  elseif weights == :cotangent
    v = vertices(ℳ)
    @assert eltype(ℳ) <: Triangle "cotangent weights only defined for triangle meshes"
    cotangentlaplacian!(L, 𝒩, v)
  else
    throw(ArgumentError("invalid discretization weights"))
  end

  L
end

function uniformlaplacian!(L, 𝒩)
  n = size(L, 1)
  for i in 1:n
    js = 𝒩(i)
    for j in js
      L[i,j] = 1 / length(js)
    end
    L[i,i] = -1
  end
end

function cotangentlaplacian!(L, 𝒩, v)
  n = size(L, 1)
  for i in 1:n
    js = CircularVector(𝒩(i))
    for k in 1:length(js)
      j₋, j, j₊ = js[k-1], js[k], js[k+1]
      vᵢ, vⱼ =  v[i],  v[j]
      v₋, v₊ = v[j₋], v[j₊]
      αᵢⱼ = ∠(vⱼ, v₋, vᵢ)
      βᵢⱼ = ∠(vᵢ, v₊, vⱼ)
      L[i,j] = cot(αᵢⱼ) + cot(βᵢⱼ)
    end
    L[i,i] = -sum(L[i,js])
  end
end