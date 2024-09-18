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

  # add centroids of elements
  ∂₂₀ = Boundary{2,0}(t)
  epts = if isnothing(method.pred)
    map(1:nelements(t)) do elem
      is = ∂₂₀(elem)
      cₒ = sum(i -> to(points[i]), is) / length(is)
      withcrs(mesh, cₒ)
    end
  else
    pts = eltype(points)[]
    for elem in mesh
      e = element(mesh, elem)
      if method.pred(e)
        push!(pts, centroid(e))
      end
    end
    pts
  end

  # original vertices
  vpts = points

  # new points in refined mesh
  newpoints = [vpts; epts]

  offset = length(vpts)

  # connect vertices into new triangles
  newconnec = Connectivity{Triangle,3}[]
  if isnothing(method.pred)
    for elem in 1:nelements(t)
      verts = ∂₂₀(elem)
      nv = length(verts)
      for i in 1:nv
        u = elem + offset
        v = verts[mod1(i, nv)]
        w = verts[mod1(i + 1, nv)]
        tri = connect((u, v, w))
        push!(newconnec, tri)
      end
    end
  else
    u = offset
    for elem in 1:nelements(t)
      if method.pred(element(mesh, elem))
        verts = ∂₂₀(elem)
        nv = length(verts)
        u += 1
        for i in 1:nv
          v = verts[mod1(i, nv)]
          w = verts[mod1(i + 1, nv)]
          tri = connect((u, v, w))
          push!(newconnec, tri)
        end
      else
        push!(newconnec, element(t, elem))
      end
    end
  end

  SimpleMesh(newpoints, newconnec)
end
