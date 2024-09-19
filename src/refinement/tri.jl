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
  topo = topology(mesh)

  # indices to refine
  rinds = if isnothing(method.pred)
    1:nelements(topo)
  else
    filter(i -> method.pred(mesh[i]), 1:nelements(topo))
  end

  # indices to preserve
  pinds = setdiff(1:nelements(topo), rinds)

  # add centroids of elements
  ∂₂₀ = Boundary{2,0}(topo)
  rpts = map(rinds) do elem
    is = ∂₂₀(elem)
    cₒ = sum(i -> to(points[i]), is) / length(is)
    withcrs(mesh, cₒ)
  end

  # original vertices
  vpts = points

  # new points in refined mesh
  newpoints = [vpts; rpts]

  # new connectivities in refined mesh
  newconnec = Connectivity{Triangle,3}[]

  # offset to new vertex indices
  offset = length(vpts)

  # connectivities of new triangles
  for (i, elem) in enumerate(rinds)
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

  # connectivities of preserved elements
  for elem in pinds
    push!(newconnec, element(topo, elem))
  end

  SimpleMesh(newpoints, newconnec)
end
