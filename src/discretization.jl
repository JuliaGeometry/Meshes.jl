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

discretize(box::Box, method::DiscretizationMethod) =
  discretize(boundary(box), method)

function discretize(polygon::Polygon{Dim,T}, method::DiscretizationMethod) where {Dim,T}
  # build bridges in case the polygon has holes,
  # i.e. reduce to a single outer boundary
  chain, dups = bridge(unique(polygon), width=2atol(T))

  # discretize using outer boundary
  mesh = discretize(chain, method)

  if isempty(dups)
    # nothing to be done, return mesh
    mesh
  else
    # remove duplicate vertices
    points = collect(vertices(mesh))
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

discretize(multi::Multi, method::DiscretizationMethod) =
  mapreduce(geom -> discretize(geom, method), merge, multi)

"""
    triangulate(geometry)

Discretize `geometry` of parametric dimension 2 into
triangles using an appropriate discretization method.
"""
triangulate(box::Box{2}) = discretize(box, Dehn1899())
triangulate(ngon::Ngon) = discretize(ngon, Dehn1899())
triangulate(poly::PolyArea) = discretize(poly, FIST())
triangulate(multi::Multi) = mapreduce(triangulate, merge, multi)

# ----------------
# IMPLEMENTATIONS
# ----------------

include("discretization/fist.jl")
include("discretization/dehn.jl")
