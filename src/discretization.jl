# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DiscretizationMethod

A method for discretizing geometries into meshes.
"""
abstract type DiscretizationMethod end

"""
    TriangulationMethod

A method for discretizing geometries into triangular meshes.
"""
abstract type TriangulationMethod <: DiscretizationMethod end

"""
    discretize(geometry, [method])

Discretize `geometry` with discretization `method`.

If the `method` is ommitted, a default algorithm is
used with a specific number of elements.
"""
function discretize end

"""
    BoundaryTriangulationMethod

A method for discretizing geometries into triangular meshes based on their boundary.
"""
abstract type BoundaryTriangulationMethod <: TriangulationMethod end

"""
    discretizewithin(boundary, method)

Discretize geometry within `boundary` with boundary discretization `method`.
"""
function discretizewithin end

# -----------
# DISCRETIZE
# -----------

discretize(geometry) = simplexify(geometry)

discretize(ball::Ball{𝔼{2}}) = discretize(ball, RegularDiscretization(50))

discretize(disk::Disk) = discretize(disk, RegularDiscretization(50))

discretize(sphere::Sphere{𝔼{3}}) = discretize(sphere, RegularDiscretization(50))

discretize(ellipsoid::Ellipsoid) = discretize(ellipsoid, RegularDiscretization(50))

discretize(torus::Torus) = discretize(torus, RegularDiscretization(50))

discretize(cyl::Cylinder) = discretize(cyl, RegularDiscretization(2, 50, 2))

discretize(cylsurf::CylinderSurface) = discretize(cylsurf, RegularDiscretization(50, 2))

discretize(consurf::ConeSurface) = discretize(consurf, RegularDiscretization(50, 2))

discretize(frustsurf::FrustumSurface) = discretize(frustsurf, RegularDiscretization(50, 2))

discretize(parsurf::ParaboloidSurface) = discretize(parsurf, RegularDiscretization(50))

discretize(multi::Multi) = mapreduce(discretize, merge, parent(multi))

discretize(tgeom::TransformedGeometry) = transform(tgeom)(discretize(parent(tgeom)))

discretize(mesh::Mesh) = mesh

# ----------
# FALLBACKS
# ----------

discretize(multi::Multi, method::DiscretizationMethod) =
  mapreduce(geom -> discretize(geom, method), merge, parent(multi))

discretize(tgeom::TransformedGeometry, method::DiscretizationMethod) =
  transform(tgeom)(discretize(parent(tgeom), method))

# -----------------
# BOUNDARY METHODS
# -----------------

discretize(geometry, method::BoundaryTriangulationMethod) = discretizewithin(boundary(geometry), method)

discretize(multi::Multi, method::BoundaryTriangulationMethod) =
  mapreduce(geom -> discretize(geom, method), merge, parent(multi))

function discretize(polygon::Polygon, method::BoundaryTriangulationMethod)
  # clean up polygon if necessary
  cpoly = polygon |> Repair(0) |> Repair(8)

  # handle degenerate polygons
  if nvertices(cpoly) == 1
    v = first(vertices(cpoly))
    points = [v, v, v]
    connec = [connect((1, 2, 3))]
    return SimpleMesh(points, connec)
  end

  # build bridges in case the polygon has holes,
  # i.e. reduce to a single outer boundary
  bpoly, dups = apply(Bridge(2atol(lentype(polygon))), cpoly)

  # discretize using outer boundary
  mesh = discretizewithin(boundary(bpoly), method)

  if isempty(dups)
    # nothing to be done, return mesh
    mesh
  else
    # remove duplicate vertices
    points = vertices(mesh)
    for (i, j) in dups
      points[i] = centroid(Segment(points[i], points[j]))
    end
    repeated = sort(last.(dups))
    deleteat!(points, repeated)

    # adjust connectivities
    elems = elements(topology(mesh))
    twin = Dict(reverse.(dups))
    rrep = reverse(repeated)
    einds = map(elems) do elem
      inds = indices(elem)
      [get(twin, ind, ind) for ind in inds]
    end
    for inds in einds
      for r in rrep
        for i in 1:length(inds)
          inds[i] > r && (inds[i] -= 1)
        end
      end
    end
    connec = [connect(ntuple(i -> inds[i], 3)) for inds in einds]

    # return mesh without duplicates
    SimpleMesh(points, connec)
  end
end

function discretizewithin(ring::Ring, method::BoundaryTriangulationMethod)
  # collect vertices to get rid of static containers
  points = collect(vertices(ring))

  # discretize within 2D ring with given method
  ring2D = Ring(_proj2D(manifold(ring), points))
  mesh = discretizewithin(ring2D, method)

  # return mesh with original points
  SimpleMesh(points, topology(mesh))
end

_proj2D(::Type{𝔼{3}}, points) = proj2D(points)

function _proj2D(::Type{🌐}, points)
  map(points) do p
    latlon = convert(LatLon, coords(p))
    flat(Point(latlon))
  end
end

# -----------
# SIMPLEXIFY
# -----------

"""
    simplexify(object)

Discretize `object` into simplices using an
appropriate discretization method.

### Notes

This function is sometimes called "triangulate"
when the `object` has parametric dimension 2.
"""
function simplexify end

simplexify(geometry) = simplexify(discretize(geometry))

simplexify(box::Box{𝔼{1}}) = SimpleMesh(collect(extrema(box)), GridTopology(1))

simplexify(box::Box{𝔼{2}}) = discretize(box, FanTriangulation())

simplexify(box::Box{𝔼{3}}) = discretize(box, ManualDiscretization())

simplexify(seg::Segment) = SimpleMesh(pointify(seg), GridTopology(1))

function simplexify(chain::Chain)
  np = nvertices(chain) + isclosed(chain)
  ip = isperiodic(chain)

  points = collect(vertices(chain))
  topo = GridTopology((np - 1,), ip)

  SimpleMesh(points, topo)
end

simplexify(bezier::BezierCurve) = discretize(bezier, RegularDiscretization(50))

simplexify(sphere::Sphere{𝔼{2}}) = discretize(sphere, RegularDiscretization(50))

simplexify(circle::Circle) = discretize(circle, RegularDiscretization(50))

simplexify(poly::Polygon) = discretize(poly, nvertices(poly) > 5000 ? DelaunayTriangulation() : DehnTriangulation())

simplexify(poly::Polyhedron) = discretize(poly, ManualDiscretization())

simplexify(multi::Multi) = mapreduce(simplexify, merge, parent(multi))

function simplexify(mesh::Mesh)
  # retrieve vertices and topology
  points = vertices(mesh)
  topo = topology(mesh)

  # check if there is something to do
  all(issimplex, elements(topo)) && return mesh

  # initialize vector of global indices
  ginds = Vector{Int}[]

  # simplexify each element and append global indices
  for connec in elements(topo)
    # materialize element and indices
    elem = materialize(connec, points)
    inds = indices(connec)

    # simplexify element
    mesh′ = simplexify(elem)
    topo′ = topology(mesh′)
    connecs′ = elements(topo′)

    # convert from local to global indices
    einds = [[inds[i] for i in indices(c′)] for c′ in connecs′]

    # save global indices
    append!(ginds, einds)
  end

  # simplex type for parametric dimension
  PL = paramdim(mesh) == 2 ? Triangle : Tetrahedron
  NV = nvertices(PL)

  # new connectivities
  newconnecs = [connect(ntuple(i -> inds[i], NV), PL) for inds in ginds]

  SimpleMesh(points, newconnecs)
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("discretization/fan.jl")
include("discretization/dehn.jl")
include("discretization/held.jl")
include("discretization/delaunay.jl")
include("discretization/manual.jl")
include("discretization/regular.jl")
