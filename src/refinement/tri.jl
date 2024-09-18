# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TriRefinement([pred])

Refinement of polygonal meshes into triangles.
A n-gon for which the predicate `pred` holds true
is subdivided into n triangles. The method refines all
n-gons if the `pred` is ommited.
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

  # indices to refine
  einds = if isnothing(method.pred)
    1:nelements(t)
  else
    einds = filter(i -> method.pred(mesh[i]), 1:nelements(t))
  end

  # indices to preserve
  vinds = setdiff(einds, 1:nelements(t))

  # add centroids of elements
  ∂₂₀ = Boundary{2,0}(t)
  epts = map(einds) do elem
    is = ∂₂₀(elem)
    cₒ = sum(i -> to(points[i]), is) / length(is)
    withcrs(mesh, cₒ)
  end

  # original vertices
  vpts = points

  # new points in refined mesh
  newpoints = [vpts; epts]

  # offset to new vertex indices
  offset = length(vpts)

  # new connectivities in refined mesh
  newconnec = Connectivity{Triangle,3}[]
  
  # push connects of preserved elements
  for elem in vinds
    push!(newconnec, element(t, elem))
  end
  
  # connect vertices into new triangles
  for (i, elem) in enumerate(einds)
    verts = ∂₂₀(elem)
    nv = length(verts)
    for j in 1:nv
      u = i + offset
      v = verts[mod1(j, nv)]
      w = verts[mod1(j + 1, nv)]
      tri = connect((u, v, w))
      push!(newconnec, tri)
    end
  end

  SimpleMesh(newpoints, newconnec)
end
