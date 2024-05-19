# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    laplacematrix(mesh; weights=:cotangent)

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
function laplacematrix(mesh; weights=:cotangent)
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
      L[i, j] = 1 / length(js)
    end
    L[i, i] = -1
  end
end

function cotangentlaplacian!(L, 𝒩, v)
  n = size(L, 1)
  for i in 1:n
    js = CircularVector(𝒩(i))
    for k in 1:length(js)
      j₋, j, j₊ = js[k - 1], js[k], js[k + 1]
      vᵢ, vⱼ = v[i], v[j]
      v₋, v₊ = v[j₋], v[j₊]
      αᵢⱼ = ∠(vⱼ, v₋, vᵢ)
      βᵢⱼ = ∠(vᵢ, v₊, vⱼ)
      L[i, j] = cot(αᵢⱼ) + cot(βᵢⱼ)
    end
    L[i, i] = -sum(L[i, js])
  end
end

"""
    measurematrix(mesh)

The measure (or "mass") matrix of the `mesh`, i.e. a diagonal
matrix with entries `Mᵢᵢ = 2Aᵢ` where `Aᵢ` is (one-third of) the
sum of the areas of triangles sharing vertex `i`.

The discrete cotangent Laplace-Beltrami operator can be written
as `Δ = M⁻¹L`. When solving systems of the form `Δu = f`, it
is useful to write `Lu = Mf` and exploit the symmetry of `L`.
"""
function measurematrix(mesh)
  # convert to half-edge topology
  ℳ = topoconvert(HalfEdgeTopology, mesh)

  # retrieve coboundary relation
  ∂ = Coboundary{0,2}(topology(ℳ))

  # pre-compute all measures
  A = measure.(ℳ)

  # initialize matrix
  n = nvertices(ℳ)
  M = oneunit(eltype(A)) * I(n)

  # fill matrix with measures
  for i in 1:n
    Aᵢ = sum(A[∂(i)]) / 3
    M[i, i] = 2Aᵢ
  end

  M
end

"""
    adjacencymatrix(mesh)

Return the adjacency matrix of the elements of the `mesh`
using the adjacency relation of the underlying topology.
"""
function adjacencymatrix(mesh)
  t = topology(mesh)
  D = paramdim(mesh)
  𝒜 = Adjacency{D}(t)

  # initialize matrix
  n = nelements(mesh)
  A = spzeros(Int, n, n)

  # fill in matrix
  for i in 1:n, j in 𝒜(i)
    A[i, j] = 1
  end

  A
end
