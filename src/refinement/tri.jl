# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TriRefinement()

Refinement of polygonal meshes into triangles.
A n-gon is subdivided into n-2 triangles.
"""
struct TriRefinement <: RefinementMethod end

function refine(mesh, ::TriRefinement)
  points = vertices(mesh)
  elems  = elements(mesh)
  topo   = topology(mesh)
  connec = elements(topo)

  # initialize vector of global indices
  ginds = Vector{Int}[]

  # triangulate each element and append global indices
  for (e, c) in zip(elems, connec)
    # triangulate single element
    mesh′   = triangulate(e)
    topo′   = topology(mesh′)
    connec′ = elements(topo′)

    # global indices
    inds = indices(c)

    # convert from local to global indices
    einds = [[inds[i] for i in indices(c′)] for c′ in connec′]

    # save global indices
    append!(ginds, einds)
  end

  # new points and connectivities
  newpoints = collect(points)
  newconnec = connect.(Tuple.(ginds), Triangle)

  SimpleMesh(newpoints, newconnec)
end
