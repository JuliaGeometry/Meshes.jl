# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Dehn1899()

Max Dehns' triangulation proved in 1899.

The algorithm is described in the first chapter of Devadoss & Rourke 2011,
and is based on a theorem derived in 1899 by the German mathematician Max Dehn.
See <https://en.wikipedia.org/wiki/Two_ears_theorem>.

Because the algorithm relies on recursion, it is mostly appropriate for polygons
with small number of vertices.

## References

* Devadoss, S & Rourke, J. 2011. [Discrete and computational geometry]
  (https://press.princeton.edu/books/hardcover/9780691145532/discrete-and-computational-geometry)
"""
struct Dehn1899 <: BoundaryDiscretizationMethod end

function discretizewithin(ring::Ring{2}, ::Dehn1899)
  # points on resulting mesh
  points = collect(vertices(ring))

  # Dehn's recursion
  connec = dehn1899(points, 1:length(points))

  SimpleMesh(points, connec)
end

function dehn1899(v::AbstractVector{<:Point}, inds)
  I = CircularVector(inds)
  n = length(I)

  if n > 3 # split chain
    # find lowerleft vertex
    i = first(sortperm(to.(v[I])))

    # left/right chains
    linds = (i - 1):(i + 1)
    rinds = (i + 1):(i + n - 1)

    # check if candidate diagonal is valid
    Δ = Triangle(v[I[i - 1]], v[I[i]], v[I[i + 1]])
    intriangle = findall(j -> v[I[j]] ∈ Δ, rinds[2:(end - 1)])
    isdiag = isempty(intriangle)

    # adjust diagonal if necessary
    if !isdiag
      l = Line(v[I[i - 1]], v[I[i + 1]])
      js = rinds[intriangle .+ 1]
      k = argmax([evaluate(Euclidean(), l, v[I[j]]) for j in js])
      j = js[k]
      linds = i:j
      rinds = j:(i + n)
    end

    # we adjust the circular indices and
    # use `inds` instead of `I` in the
    # recursion to avoid memory copies
    linds = [mod1(ind, n) for ind in linds]
    rinds = [mod1(ind, n) for ind in rinds]

    # perform recursion
    left = dehn1899(v, inds[linds])
    right = dehn1899(v, inds[rinds])
    [left; right]
  else
    # return the triangle
    [connect(Tuple(inds), Triangle)]
  end
end
