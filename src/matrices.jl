# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# helper function to select default Laplacian discretization
laplacekind(mesh) = eltype(mesh) <: Triangle ? :cotangent : :uniform

# helper function to convert topology if necessary
adjusttopo(topo::SimpleTopology) = convert(HalfEdgeTopology, topo)
adjusttopo(topo) = topo

"""
    laplacematrix(mesh; kind=nothing)

The Laplace-Beltrami (a.k.a. Laplacian) matrix of the `mesh`.
Optionally, specify the `kind` of discretization.

## Available discretizations

* `:uniform`   - `Lᵢⱼ = 1 / |𝒜(i)|, ∀j ∈ 𝒜(i)`
* `:cotangent` - `Lᵢⱼ = cot(αᵢⱼ) + cot(βᵢⱼ), ∀j ∈ 𝒜(i)`

where `𝒜(i)` is the adjacency relation at vertex `i`.

## References

* Botsch et al. 2010. [Polygon Mesh Processing](http://www.pmp-book.org).

* Pinkall, U. & Polthier, K. 1993. [Computing discrete minimal surfaces and their conjugates]
  (https://projecteuclid.org/journals/experimental-mathematics/volume-2/issue-1/Computing-discrete-minimal-surfaces-and-their-conjugates/em/1062620735.full).
"""
function laplacematrix(mesh; kind=nothing)
  # select default discretization
  𝒦 = isnothing(kind) ? laplacekind(mesh) : kind

  # sanity checks
  𝒦 == :cotangent && assertion(eltype(mesh) <: Triangle, "cotangent weights only defined for triangle meshes")

  # adjust topology if necessary
  𝒯 = adjusttopo(topology(mesh))

  # retrieve adjacency relation
  𝒜 = Adjacency{0}(𝒯)

  # initialize matrix
  n = nvertices(mesh)
  L = spzeros(n, n)

  # fill matrix
  if 𝒦 == :uniform
    uniformlaplacian!(L, 𝒜)
  elseif 𝒦 == :cotangent
    cotangentlaplacian!(L, 𝒜, vertices(mesh))
  end

  L
end

function uniformlaplacian!(L, 𝒜)
  n = size(L, 1)
  for i in 1:n
    js = 𝒜(i)
    for j in js
      L[i, j] = 1 / length(js)
    end
    L[i, i] = -1
  end
end

function cotangentlaplacian!(L, 𝒜, v)
  n = size(L, 1)
  for i in 1:n
    js = 𝒜(i)
    m = length(js)
    for k in 1:m
      j₋ = js[mod1(k - 1, m)]
      j = js[mod1(k, m)]
      j₊ = js[mod1(k + 1, m)]
      vᵢ, vⱼ = v[i], v[j]
      v₋, v₊ = v[j₋], v[j₊]
      αᵢⱼ = ∠(vⱼ, v₋, vᵢ)
      βᵢⱼ = ∠(vᵢ, v₊, vⱼ)
      L[i, j] = cot(αᵢⱼ) + cot(βᵢⱼ)
    end
    L[i, i] = -sum(j -> L[i, j], js)
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
  # adjust topology if necessary
  𝒯 = adjusttopo(topology(mesh))

  # parametric dimension
  D = paramdim(mesh)

  # retrieve coboundary relation
  𝒞 = Coboundary{0,D}(𝒯)

  # pre-compute all measures
  A = map(measure, mesh)

  # initialize matrix
  n = nvertices(mesh)
  M = oneunit(eltype(A)) * I(n)

  # fill matrix
  for i in 1:n
    js = 𝒞(i)
    Aᵢ = sum(j -> A[j], js) / 3
    M[i, i] = 2Aᵢ
  end

  M
end

"""
    adjacencymatrix(mesh; rank=paramdim(mesh))

The adjacency matrix of the `mesh` using the adjacency
relation of given `rank` for the underlying topology.
"""
function adjacencymatrix(mesh; rank=paramdim(mesh))
  # adjust topology if necessary
  𝒯 = adjusttopo(topology(mesh))

  # retrieve adjacency relation
  𝒜 = Adjacency{rank}(𝒯)

  # initialize matrix
  n = nfaces(mesh, rank)
  A = spzeros(Int, n, n)

  # fill in matrix
  for i in 1:n, j in 𝒜(i)
    A[i, j] = 1
  end

  A
end
