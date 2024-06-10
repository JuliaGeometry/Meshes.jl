# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DiscretizationMethod

A method for discretizing geometries into meshes.
"""
abstract type DiscretizationMethod end

"""
    discretize(geometry, [method])

Discretize `geometry` with discretization `method`.

If the `method` is ommitted, a default algorithm is
used with a specific number of elements.
"""
function discretize end

"""
    BoundaryDiscretizationMethod

A method for discretizing geometries based on their boundary.
"""
abstract type BoundaryDiscretizationMethod <: DiscretizationMethod end

"""
    discretizewithin(boundary, method)

Discretize geometry within `boundary` with boundary discretization `method`.
"""
function discretizewithin end

discretize(geometry::Geometry, method::BoundaryDiscretizationMethod) = discretizewithin(boundary(geometry), method)

function discretize(polygon::Polygon, method::BoundaryDiscretizationMethod)
  # clean up polygon if necessary
  cpoly = polygon |> Repair{0}() |> Repair{8}()

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
    connec = connect.(Tuple.(einds))

    # return mesh without duplicates
    SimpleMesh(points, connec)
  end
end

discretize(multi::Multi, method::BoundaryDiscretizationMethod) =
  mapreduce(geom -> discretize(geom, method), merge, parent(multi))

function discretizewithin(ring::Ring{3}, method::BoundaryDiscretizationMethod)
  # collect vertices to get rid of static containers
  points = collect(vertices(ring))

  # discretize within 2D ring with given method
  ring2D = Ring(proj2D(points))
  mesh = discretizewithin(ring2D, method)

  # return mesh with original points
  SimpleMesh(points, topology(mesh))
end

# ----------------
# DEFAULT METHODS
# ----------------

discretize(geometry) = simplexify(geometry)

discretize(ball::Ball{2}) = discretize(ball, RegularDiscretization(50))

discretize(disk::Disk) = discretize(disk, RegularDiscretization(50))

discretize(sphere::Sphere{3}) = discretize(sphere, RegularDiscretization(50))

discretize(ellipsoid::Ellipsoid) = discretize(ellipsoid, RegularDiscretization(50))

discretize(torus::Torus) = discretize(torus, RegularDiscretization(50))

discretize(cylsurf::CylinderSurface) = discretize(cylsurf, RegularDiscretization(50, 2))

discretize(consurf::ConeSurface) = discretize(consurf, RegularDiscretization(50, 2))

discretize(frustsurf::FrustumSurface) = discretize(frustsurf, RegularDiscretization(50, 2))

discretize(parsurf::ParaboloidSurface) = discretize(parsurf, RegularDiscretization(50))

discretize(multi::Multi) = mapreduce(discretize, merge, parent(multi))

discretize(mesh::Mesh) = mesh

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

simplexify(box::Box{1}) = SimpleMesh(collect(extrema(box)), GridTopology(1))

simplexify(seg::Segment) = SimpleMesh(pointify(seg), GridTopology(1))

function simplexify(chain::Chain)
  np = nvertices(chain) + isclosed(chain)
  ip = isperiodic(chain)

  points = collect(vertices(chain))
  topo = GridTopology((np - 1,), ip)

  SimpleMesh(points, topo)
end

simplexify(bezier::BezierCurve) = discretize(bezier, RegularDiscretization(50))

simplexify(sphere::Sphere{2}) = discretize(sphere, RegularDiscretization(50))

simplexify(circle::Circle) = discretize(circle, RegularDiscretization(50))

simplexify(box::Box{2}) = discretize(box, FanTriangulation())

simplexify(box::Box{3}) = discretize(box, Tetrahedralization())

simplexify(poly::Polygon) = discretize(poly, nvertices(poly) > 5000 ? HeldTriangulation() : DehnTriangulation())

simplexify(poly::Polyhedron) = discretize(poly, Tetrahedralization())

simplexify(multi::Multi) = mapreduce(simplexify, merge, parent(multi))

function simplexify(mesh::Mesh)
  points = vertices(mesh)
  elems = elements(mesh)
  topo = topology(mesh)
  connec = elements(topo)

  # initialize vector of global indices
  ginds = Vector{Int}[]

  # simplexify each element and append global indices
  for (e, c) in zip(elems, connec)
    # simplexify single element
    mesh′ = simplexify(e)
    topo′ = topology(mesh′)
    connec′ = elements(topo′)

    # global indices
    inds = indices(c)

    # convert from local to global indices
    einds = [[inds[i] for i in indices(c′)] for c′ in connec′]

    # save global indices
    append!(ginds, einds)
  end

  # simplex type for parametric dimension
  PL = paramdim(mesh) == 2 ? Triangle : Tetrahedron

  # new connectivities
  newconnec = connect.(Tuple.(ginds), PL)

  SimpleMesh(points, newconnec)
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("discretization/fan.jl")
include("discretization/dehn.jl")
include("discretization/held.jl")
include("discretization/tetra.jl")
include("discretization/regular.jl")
