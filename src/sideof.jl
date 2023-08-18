# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    IntersectionType

The different types of sides that a point may lie in relation to a
geometry or domain. Type `SideType` in a Julia session to see the
full list.
"""
@enum SideType begin
  IN
  OUT
  ON
  LEFT
  RIGHT
end

"""
    sideof(point, line)

Determines on which side the `point` is in relation to the `line`.
Possible results are `LEFT`, `RIGHT` or `ON` the `line`.

### Notes

* Assumes the orientation of `Segment(line(0), line(1))`.
"""
function sideof(p::Point{2,T}, l::Line{2,T}) where {T}
  area = signarea(p, l(0), l(1))
  ifelse(area > atol(T), LEFT, ifelse(area < -atol(T), RIGHT, ON))
end

"""
    sideof(point, ring)

Determines on which side of the `ring` the `point` lies.
Possible results are `IN` or `OUT` the `ring`.
"""
function sideof(p::Point{2,T}, r::Ring{2,T}) where {T}
  w = windingnumber(p, r)
  ifelse(isapprox(w, zero(T), atol=atol(T)), OUT, IN)
end

"""
    sideof(point, mesh)

Determines on which side of the surface `mesh` the `point` lies.
Possible results are `IN`, `OUT`, or `ON` the `mesh`.

### Notes

* Uses a ray-casting algorithm.
"""
function sideof(point::Point{3,T}, mesh::Mesh{3,T}) where {T}
  @assert paramdim(mesh) == 2 "sideof only defined for surface meshes"
  (eltype(mesh) <: Triangle) || return sideof(point, simplexify(mesh))

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
