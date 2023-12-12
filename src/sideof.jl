# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    IntersectionType

The different types of sides that a point may lie in relation to a
boundary geometry or mesh. Type `SideType` in a Julia session to see
the full list.
"""
@enum SideType begin
  IN
  OUT
  ON
  LEFT
  RIGHT
end

"""
    sideof(points, object)

Determines on which side the `points` are in relation to the geometric
`object`, which can be a boundary `geometry` or `mesh`.
"""
function sideof end

# ---------
# GEOMETRY
# ---------

"""
    sideof(point, line)

Determines on which side the `point` is in relation to the `line`.
Possible results are `LEFT`, `RIGHT` or `ON` the `line`.

### Notes

* Assumes the orientation of `Segment(line(0), line(1))`.
"""
function sideof(point::Point{2,T}, line::Line{2,T}) where {T}
  a = signarea(point, line(0), line(1))
  ifelse(a > atol(T), LEFT, ifelse(a < -atol(T), RIGHT, ON))
end

"""
    sideof(point, ring)

Determines on which side the `point` is in relation to the `ring`.
Possible results are `IN` or `OUT` the `ring`.
"""
function sideof(point::Point{2,T}, ring::Ring{2,T}) where {T}
  w = windingnumber(point, ring)
  ifelse(isapprox(w, zero(T), atol=atol(T)), OUT, IN)
end

# fallback for iterable of points
sideof(points, geom::Geometry) = map(point -> sideof(point, geom), points)

# -----
# MESH
# -----

"""
    sideof(point, mesh)

Determines on which side the `point` is in relation to the surface `mesh`.
Possible results are `IN`, `OUT`, or `ON` the `mesh`.
"""
sideof(point::Point{3}, mesh::Mesh{3}) = sideof([point], mesh) |> first

function sideof(points, mesh::Mesh{3})
  @assert paramdim(mesh) == 2 "sideof only defined for surface meshes"

  # triangulate mesh if necessary
  (eltype(mesh) <: Triangle) || return sideof(points, simplexify(mesh))

  # SVD basis for projection
  verts = vertices(mesh)
  basis = svdbasis(verts)

  # project points and mesh on SVD basis
  ppts = proj(points, basis)
  pmsh = SimpleMesh(proj(verts, basis), topology(mesh))

  # bounding box of projected triangles
  boxs = boundingbox.(pmsh)

  # retreive basis
  u, v = basis

  # loop over query points
  map(1:length(points)) do i
    # find triangles with intersecting bounding box
    inds = findall(b -> ppts[i] ∈ b, boxs)

    if isempty(inds)
      OUT # the point is outside the mesh
    else
      # filter original triangle mesh
      fmsh = view(mesh, inds)

      # define ray from SVD basis
      nray = Ray(points[i], u × v)

      # ray casting algorithm with filtered mesh
      _raycasting(nray, fmsh)
    end
  end
end

function _raycasting(ray::Ray{3,T}, tris) where {T}
  ecross = 0 # number of edge crosses
  inters = false # do we have intersection?
  points = Point{3,T}[]
  for tri in tris
    I = intersection(ray, tri)
    if type(I) == Crossing
      inters = !inters
    elseif type(I) ∈ (EdgeTouching, CornerTouching, Touching)
      return ON
    elseif type(I) == EdgeCrossing
      ecross += 1
    elseif type(I) == CornerCrossing
      point = get(I)
      if !any(≈(point), points)
        push!(points, point)
        inters = !inters
      end
    end
  end

  # check how many edges we crossed
  isodd(ecross ÷ 2) && (inters = !inters)

  inters ? IN : OUT
end
