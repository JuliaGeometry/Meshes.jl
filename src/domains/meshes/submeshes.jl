# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SubMesh{M,CRS}

A submesh of geometries in a given manifold `M` with point
coordinates specified in a coordinate reference system `CRS`.
"""
const SubMesh{M<:Manifold,C<:CRS} = SubDomain{M,C,<:Mesh{M,C}}

function materialize(m::SubMesh)
  # retrieve parent mesh and indices
  pmesh = parent(m)
  pinds = parentindices(m)

  # indices of other elements
  oinds = setdiff(1:nelements(pmesh), pinds)

  # topology and vertices of parent mesh
  ptopo = topology(pmesh)
  pvert = vertices(pmesh)

  # find indices of vertices that are not in view
  select = [i for pind in pinds for i in indices(element(ptopo, pind))]
  notsel = [i for oind in oinds for i in indices(element(ptopo, oind))]
  others = setdiff(notsel, select)

  # map old vertex indices to new vertex indices
  newind = collect(1:nvertices(ptopo))
  @inbounds for i in others
    newind[i+1:end] .-= 1
  end
  newind[others] .= 0

  # map new vertex indices to old vertex indices
  indmax = 0
  oldind = Dict{Int,Int}()
  for (i, n) in enumerate(newind)
    n > indmax && (indmax = n)
    n > 0 && (oldind[n] = i)
  end

  # retrieve new vertices from parent vertices
  points = [pvert[oldind[n]] for n in 1:indmax]

  # construct new connectivities with updated indices
  connec = @inbounds map(pinds) do pind
    elem = element(ptopo, pind)
    inds = map(i -> newind[i], indices(elem))
    connect(inds, pltype(elem))
  end

  SimpleMesh(points, connec)
end

"""
    SubGrid{M,CRS}

A subgrid of geometries in a given manifold `M` with point
coordinates specified in a coordinate reference system `CRS`.
"""
const SubGrid{M<:Manifold,C<:CRS} = SubDomain{M,C,<:Grid{M,C}}
