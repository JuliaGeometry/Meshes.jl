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

function discretize(geometry::TransformedGeometry)
  T = numtype(lentype(geometry))
  mesh = if hasdistortedboundary(geometry)
    discretize(parent(geometry), MaxLengthDiscretization(T(1000) * u"km"))
  else
    discretize(parent(geometry))
  end
  transform(geometry)(mesh)
end

discretize(mesh::Mesh) = mesh

# ----------
# FALLBACKS
# ----------

discretize(multi::Multi, method::DiscretizationMethod) =
  mapreduce(geom -> discretize(geom, method), merge, parent(multi))

discretize(geometry::TransformedGeometry, method::DiscretizationMethod) =
  transform(geometry)(discretize(parent(geometry), method))

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
  # retrieve vertices of ring
  points = collect(eachvertex(ring))

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
include("discretization/adaptive.jl")
