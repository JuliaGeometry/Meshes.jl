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

If the `method` is ommitted, a default is used as a
function of the `geometry`. Geometries over the `🌐`
manifold are refined until the segments are shorter
than a maximum length.
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

discretize(geometry::Geometry) = simplexify(geometry)

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

function discretize(geometry::TransformedGeometry)
  pmesh = discretize(parent(geometry))
  tmesh = pmesh |> transform(geometry)
  _mayberefinemaxlen(pmesh, tmesh)
end

discretize(mesh::Mesh) = mesh

# ----------
# FALLBACKS
# ----------

function discretize(geometry::Geometry, method::DiscretizationMethod)
  if manifold(geometry) == 🌐
    _discretize🌐(geometry, method)
  else
    _discretize𝔼(geometry, method)
  end
end

_discretize🌐(geometry::Geometry, method::DiscretizationMethod) = _refinemaxlen(_discretize(geometry, method))

_discretize𝔼(geometry::Geometry, method::DiscretizationMethod) = _discretize(geometry, method)

discretize(multi::Multi, method::DiscretizationMethod) =
  mapreduce(geom -> discretize(geom, method), merge, parent(multi))

function discretize(geometry::TransformedGeometry, method::DiscretizationMethod)
  pmesh = discretize(parent(geometry), method)
  tmesh = pmesh |> transform(geometry)
  _mayberefinemaxlen(pmesh, tmesh)
end

# -----------------
# BOUNDARY METHODS
# -----------------

discretize(geometry::Geometry, method::BoundaryTriangulationMethod) = discretizewithin(boundary(geometry), method)

function discretize(polygon::Polygon, method::BoundaryTriangulationMethod)
  # clean up polygon if necessary
  cpoly = polygon |> Repair(0) |> Repair(8)

  # handle degenerate polygons
  if nvertices(cpoly) == 1
    points = fill(vertex(cpoly, 1), 3)
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
  if manifold(ring) == 🌐
    _discretizewithin🌐(ring, method)
  else
    _discretizewithin𝔼(ring, method)
  end
end

function _discretizewithin🌐(ring::Ring, method::BoundaryTriangulationMethod)
  points = collect(eachvertex(ring))
  ring𝔼2 = Ring([flat(p |> Proj(LatLon)) for p in points])
  mesh𝔼2 = _discretizewithin𝔼2(ring𝔼2, method)
  mesh🌐 = SimpleMesh(points, topology(mesh𝔼2))
  _refinemaxlen(mesh🌐)
end

function _discretizewithin𝔼(ring::Ring, method::BoundaryTriangulationMethod)
  if embeddim(ring) == 2
    _discretizewithin𝔼2(ring, method)
  elseif embeddim(ring) == 3
    points = collect(eachvertex(ring))
    ring𝔼2 = Ring(proj2D(points))
    mesh𝔼2 = _discretizewithin𝔼2(ring𝔼2, method)
    SimpleMesh(points, topology(mesh𝔼2))
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

simplexify(geometry::Geometry) = simplexify(discretize(geometry))

simplexify(box::Box) = discretize(box, ManualSimplexification())

simplexify(chain::Chain) = discretize(chain, ManualSimplexification())

simplexify(bezier::BezierCurve) = discretize(bezier, RegularDiscretization(50))

simplexify(curve::ParametrizedCurve) = discretize(curve, RegularDiscretization(50))

simplexify(sphere::Sphere{𝔼{2}}) = discretize(sphere, RegularDiscretization(50))

simplexify(circle::Circle) = discretize(circle, RegularDiscretization(50))

simplexify(tri::Triangle) = discretize(tri, ManualSimplexification())

simplexify(poly::Polygon) = discretize(poly, nvertices(poly) > 5000 ? DelaunayTriangulation() : DehnTriangulation())

simplexify(poly::Polyhedron) = discretize(poly, ManualSimplexification())

simplexify(multi::Multi) = mapreduce(simplexify, merge, parent(multi))

function simplexify(geometry::TransformedGeometry)
  pmesh = simplexify(parent(geometry))
  tmesh = pmesh |> transform(geometry)
  _mayberefinemaxlen(pmesh, tmesh)
end

function simplexify(mesh::Mesh)
  # retrieve vertices and connectivities
  points = vertices(mesh)
  connec = elements(topology(mesh))

  # check if there is something to do
  all(issimplex, connec) && return mesh

  # function barrier for optimal performance
  _simplexify(points, connec)
end

function _simplexify(points, connec)
  # initialize vector of global indices
  ginds = Vector{Int}[]

  # simplexify each element and append global indices
  for c in connec
    # manually union-split most common polytopes
    # for type stability and maximum performance
    if c isa Connectivity{Triangle,3}
      _appendinds!(ginds, c, points)
    elseif c isa Connectivity{Quadrangle,4}
      _appendinds!(ginds, c, points)
    elseif c isa Connectivity{Tetrahedron,4}
      _appendinds!(ginds, c, points)
    elseif c isa Connectivity{Hexahedron,8}
      _appendinds!(ginds, c, points)
    else
      _appendinds!(ginds, c, points)
    end
  end

  # new connectivities
  newconnec = _newconnec(ginds, connec)

  SimpleMesh(points, newconnec)
end

function _appendinds!(ginds, connec, points)
  # materialize element and indices
  elem = materialize(connec, points)
  inds = indices(connec)

  # simplexify element
  mesh = simplexify(elem)
  topo = topology(mesh)

  # convert from local to global indices
  einds = [[inds[i] for i in indices(c)] for c in elements(topo)]

  # save global indices
  append!(ginds, einds)
end

_newconnec(ginds, connec) = _newconnec(ginds, Val(paramdim(first(connec))))
_newconnec(ginds, ::Val{2}) = [connect(ntuple(i -> inds[i], 3), Triangle) for inds in ginds]
_newconnec(ginds, ::Val{3}) = [connect(ntuple(i -> inds[i], 4), Tetrahedron) for inds in ginds]

# ----------------
# IMPLEMENTATIONS
# ----------------

include("discretization/fan.jl")
include("discretization/dehn.jl")
include("discretization/held.jl")
include("discretization/delaunay.jl")
include("discretization/manual.jl")
include("discretization/regular.jl")
include("discretization/maxlength.jl")

# -----------------
# HELPER FUNCTIONS
# -----------------

function _mayberefinemaxlen(pmesh, tmesh)
  # if the manifold changes from 🌐 to 𝔼 or vice-versa
  # the mesh might be distorted in the target manifold
  # we refine the mesh further until the segments have
  # a maximum predefined length in physical units
  Mₚ, Mₜ = manifold(pmesh), manifold(tmesh)
  changed = (Mₚ == 🌐 && Mₜ != 🌐) || (Mₚ != 🌐 && Mₜ == 🌐)
  changed ? _refinemaxlen(tmesh) : tmesh
end

function _refinemaxlen(tmesh)
  T = numtype(lentype(tmesh))
  l = T(1000) * u"km"
  refine(tmesh, MaxLengthRefinement(l))
end
