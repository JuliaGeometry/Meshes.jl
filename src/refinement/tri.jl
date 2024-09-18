# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TriRefinement()

Refinement of polygonal meshes into triangles.
A n-gon is subdivided into n triangles.

    TriRefinement(pred)

Refine only the elements of the mesh where
`pred(element)` returns `true`.
"""
struct TriRefinement{F} <: RefinementMethod
  pred::F
end

TriRefinement() = TriRefinement(nothing)

function refine(mesh, method::TriRefinement)
  assertion(paramdim(mesh) == 2, "TriRefinement only defined for surface meshes")
  (eltype(mesh) <: Triangle) || return simplexify(mesh)

  # retrieve geometry and topology
  points = vertices(mesh)
  connec = topology(mesh)

  # convert to half-edge structure
  t = convert(HalfEdgeTopology, connec)

  # used to extract the vertex indices
  ∂₂₀ = Boundary{2,0}(t)

  # offset to new vertex indices
  offset = length(points)

  # predicate function
  pred = if isnothing(method.pred)
    _ -> true
  else
    eᵢ -> method.pred(element(mesh, eᵢ))
  end

  # add centroids of elements and connect vertices 
  # into new triangles if necessary
  newpoints = copy(points)
  newconnec = Connectivity{Triangle,3}[]
  for eᵢ in 1:nelements(t)
    # check if the element should be refined
    if pred(eᵢ)
      verts = ∂₂₀(eᵢ)
      nv = length(verts)
  
      # add new centroid vertex
      cₒ = sum(i -> to(points[i]), verts) / length(verts)
      pₒ = withcrs(mesh, cₒ)
      push!(newpoints, pₒ)

      # add new connectivities
      for i in 1:nv
        u = eᵢ + offset
        v = verts[mod1(i, nv)]
        w = verts[mod1(i + 1, nv)]
        tri = connect((u, v, w))
        push!(newconnec, tri)
      end
    else
      # otherwise, just use the original connectivity
      push!(newconnec, element(t, eᵢ))
    end
  end

  SimpleMesh(newpoints, newconnec)
end
