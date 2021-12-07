# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DiscretizationMethod

A method for discretizing geometries into meshes.
"""
abstract type DiscretizationMethod end

"""
    discretize(geometry, method)

Discretize `geometry` with discretization `method`.
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

discretize(geometry, method::BoundaryDiscretizationMethod) =
  discretizewithin(boundary(geometry), method)

function discretize(polygon::Polygon{Dim,T}, method::BoundaryDiscretizationMethod) where {Dim,T}
  # build bridges in case the polygon has holes,
  # i.e. reduce to a single outer boundary
  chain, dups = bridge(unique(polygon), width=2atol(T))

  # discretize using outer boundary
  mesh = discretizewithin(chain, method)

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
    twin  = Dict(reverse.(dups))
    rrep  = reverse(repeated)
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
  mapreduce(geom -> discretize(geom, method), merge, multi)

function discretizewithin(chain::Chain{3}, method::BoundaryDiscretizationMethod)
  # collect vertices to get rid of static containers
  points = vertices(chain) |> collect

  # project points on 2D plane using SVD
  # https://math.stackexchange.com/a/99317
  X = mapreduce(coordinates, hcat, points)
  μ = sum(X, dims=2) / size(X, 2)
  Z = X .- μ
  U = svd(Z).U
  u = U[:,1]
  v = U[:,2]

  # projected points
  projected = [Point(z⋅u, z⋅v) for z in eachcol(Z)]

  # discretize within 2D chain
  chain2D = Chain([projected; first(projected)])
  mesh    = discretizewithin(chain2D, method)

  # return mesh with original points
  SimpleMesh(points, topology(mesh))
end

"""
    triangulate(object)

Triangulate `object` of parametric dimension 2 into
triangles using an appropriate triangulation method.
"""
function triangulate end

triangulate(box::Box{2}) = discretize(box, FanTriangulation())

triangulate(tri::Triangle) = discretize(tri, FanTriangulation())

triangulate(quad::Quadrangle) = discretize(quad, FanTriangulation())

triangulate(ngon::Ngon) = discretize(ngon, Dehn1899())

triangulate(poly::PolyArea) = discretize(poly, FIST())

triangulate(multi::Multi) = mapreduce(triangulate, merge, multi)

function triangulate(mesh::Mesh)
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

  # new connectivities
  newconnec = connect.(Tuple.(ginds), Triangle)

  SimpleMesh(points, newconnec)
end

triangulate(sphere::Sphere{3}) =
  discretize(sphere, RegularDiscretization(100)) |> triangulate

triangulate(ball::Ball{2}) =
  discretize(ball, RegularDiscretization(100)) |> triangulate

# ----------------
# IMPLEMENTATIONS
# ----------------

include("discretization/fan.jl")
include("discretization/regular.jl")
include("discretization/fist.jl")
include("discretization/dehn.jl")
