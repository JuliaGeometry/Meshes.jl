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
sideof(point::Point{3}, mesh::Mesh{3}) = sideof((point,), mesh) |> first

function sideof(points, mesh::Mesh{3})
  @assert paramdim(mesh) == 2 "sideof only defined for surface meshes"
  (eltype(mesh) <: Triangle) || return sideof(points, simplexify(mesh))
  map(point -> _sideof(point, mesh), points)
end

function _sideof(point::Point{3,T}, mesh::Mesh{3,T}) where {T}
  z = last.(coordinates.(extrema(mesh)))
  r = Ray(point, Vec(zero(T), zero(T), 2 * (z[2] - z[1])))

  hasintersect = false
  edgecrosses = 0
  ps = Point{3,T}[]
  for t in mesh
    I = intersection(r, t)
    if type(I) == Crossing
      hasintersect = !hasintersect
    elseif type(I) ∈ (EdgeTouching, CornerTouching, Touching)
      return ON
    elseif type(I) == EdgeCrossing
      edgecrosses += 1
    elseif type(I) == CornerCrossing
      p = get(I)
      if !any(≈(p), ps)
        push!(ps, p)
        hasintersect = !hasintersect
      end
    end
  end

  # check how many edges we crossed
  isodd(edgecrosses ÷ 2) && (hasintersect = !hasintersect)

  hasintersect ? IN : OUT
end
