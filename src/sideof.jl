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
function sideof(point::Point, line::Line)
  a = signarea(point, line(0), line(1))
  ifelse(a > atol(a), LEFT, ifelse(a < -atol(a), RIGHT, ON))
end

"""
    sideof(point, ring)

Determines on which side the `point` is in relation to the `ring`.
Possible results are `IN`, `OUT` or `ON` the `ring`.

## References

* Hao et al. 2018. [Optimal Reliable Point-in-Polygon Test and
  Differential Coding Boolean Operations on Polygons]
  (https://www.mdpi.com/2073-8994/10/10/477)
"""
function sideof(point::Point, ring::Ring)
  assertion(CoordRefSystems.ncoords(crs(point)) == 2, "points must have 2 coordinates")

  n = nvertices(ring)
  if n <= 1000 || Threads.nthreads() == 1
    _sideofserial(point, ring)
  else
    _sideofthreads(point, ring)
  end
end

function _sideofthreads(point::Point, r::Ring)
  v = vertices(r)
  k = Threads.Atomic{Int}(0)
  ison = Threads.Atomic{Bool}(false)
  Threads.@threads for i in eachindex(v)
    _ison, increment_k = _sideofcore(point, v[i], v[i + 1])
    (ison[] = _ison) && break
    increment_k && Threads.atomic_add!(k, 1)
  end
  ison[] ? ON : (iseven(k[]) ? OUT : IN)
end

function _sideofserial(point::Point, r::Ring)
  v = vertices(r)
  k = 0
  for i in eachindex(v)
    ison, increment_k = _sideofcore(point, v[i], v[i + 1])
    ison && return ON
    k += increment_k
  end
  iseven(k) ? OUT : IN
end

function _sideofcore(point::Point, pᵢ::Point, pⱼ::Point)
  # flat coordinates of query point
  p = flat(coords(point))
  xₚ, yₚ = p.x, p.y

  # Create possible return values for readability
  ISON = (true, false) # ison=true, increment_k=false
  INCREMENT_K = (false, true) # ison=false, increment_k=true
  UNCHANGED = (false, false) # ison=false, increment_k=false

  # flat coordinates of segment i -- i+1
  pᵢ = flat(coords(pᵢ))
  pⱼ = flat(coords(pⱼ))
  xᵢ, yᵢ = pᵢ.x, pᵢ.y
  xⱼ, yⱼ = pⱼ.x, pⱼ.y

  v₁ = yᵢ - yₚ
  v₂ = yⱼ - yₚ

  if (isnegative(v₁) && isnegative(v₂)) || (ispositive(v₁) && ispositive(v₂))
    # case 11, 26
    return UNCHANGED
  end

  u₁ = xᵢ - xₚ
  u₂ = xⱼ - xₚ

  if ispositive(v₂) && isnonpositive(v₁)
    # case 3, 9, 16, 21, 13, 24
    f = u₁ * v₂ - u₂ * v₁
    if ispositive(f)
      # case 3, 9
      return INCREMENT_K
    elseif isequalzero(f)
      # case 16, 21
      return ISON
    end
  elseif ispositive(v₁) && isnonpositive(v₂)
    # case 4, 10, 19, 20, 12, 25
    f = u₁ * v₂ - u₂ * v₁
    if isnegative(f)
      # case 4, 10
      return INCREMENT_K
    elseif isequalzero(f)
      # case 19, 20
      return ISON
    end
  elseif isequalzero(v₂) && isnegative(v₁)
    # case 7, 14, 17
    f = u₁ * v₂ - u₂ * v₁
    if isequalzero(f)
      # case 17
      return ISON
    end
  elseif isequalzero(v₁) && isnegative(v₂)
    # case 8, 15, 18
    f = u₁ * v₂ - u₂ * v₁
    if isequalzero(f)
      # case 18
      return ISON
    end
  elseif isequalzero(v₁) && isequalzero(v₂)
    # case 1, 2, 5, 6, 22, 23
    if isnonpositive(u₂) && isnonnegative(u₁)
      # case 1
      return ISON
    elseif isnonpositive(u₁) && isnonnegative(u₂)
      # case 2
      return ISON
    end
  end
  # case 5, 6, 7, 8, 12, 13, 14, 15, 22, 23, 24, 25
  return UNCHANGED
end

# -----
# MESH
# -----

"""
    sideof(point, mesh)

Determines on which side the `point` is in relation to the surface `mesh`.
Possible results are `IN` or `OUT` the `mesh`.
"""
sideof(point::Point, mesh::Mesh) = sideof((point,), mesh) |> first

# ----------
# FALLBACKS
# ----------

sideof(points, line::Line) = map(point -> sideof(point, line), points)

function sideof(points, object::GeometryOrDomain)
  bbox = boundingbox(object)
  isin = [point ∈ bbox for point in points]
  inds = findall(isin)
  side = fill(OUT, length(isin))
  side[inds] .= sidewithinbox(collectat(points, inds), object)
  side
end

sidewithinbox(points, ring::Ring) = map(point -> sideof(point, ring), points)

sidewithinbox(points, mesh::Mesh) = map(w -> ifelse(isapproxzero(w), OUT, IN), winding(points, mesh))
