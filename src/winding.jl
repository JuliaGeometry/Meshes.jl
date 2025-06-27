# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    winding(points, object)

Generalized winding number of `points` with respect to the geometric `object`.

## References

* Barill et al. 2018. [Fast winding numbers for soups and clouds]
  (https://dl.acm.org/doi/10.1145/3197517.3201337)
* Jacobson et al. 2013. [Robust inside-outside segmentation using generalized winding numbers]
  (https://dl.acm.org/doi/10.1145/2461912.2461916)
* Oosterom, A. & Strackee, J. 1983. [The Solid Angle of a Plane Triangle]
  (https://ieeexplore.ieee.org/document/4121581)
"""
function winding end

# ------
# RINGS
# ------

function winding(points, ring::Ring)
  v = vertices(ring)
  n = nvertices(ring)

  function w(p)
    Σ = sum(∠(flat(v[i]), flat(p), flat(v[i + 1])) for i in 1:n)
    Σ / oftype(Σ, 2π)
  end

  [w(p) for p in points]
end

winding(point::Point, ring::Ring) = winding((point,), ring) |> first

# -------
# MESHES
# -------

# Jacobson et al 2013.
function winding(points, mesh::Mesh)
  assertion(paramdim(mesh) == 2, "winding number only defined for surface meshes")
  (eltype(mesh) <: Triangle) || return winding(points, simplexify(mesh))

  function windfun(p)
    ∑ = sum(1:nelements(mesh)) do i
      v = vertices(mesh[i])
      a⃗ = v[1] - p
      b⃗ = v[2] - p
      c⃗ = v[3] - p
      a = norm(a⃗)
      b = norm(b⃗)
      c = norm(c⃗)
      n = det([a⃗ b⃗ c⃗])
      d = a * b * c + (a⃗ ⋅ b⃗) * c + (b⃗ ⋅ c⃗) * a + (c⃗ ⋅ a⃗) * b
      2atan(n, d)
    end
    ∑ / oftype(∑, 4π)
  end

  if isthreaded(points isa AbstractVector && length(points) > 10)
    _windingthread(windfun, points)
  else
    _windingserial(windfun, points)
  end
end

_windingserial(windfun, points) = map(windfun, points)

function _windingthread(windfun, points)
  T = numtype(lentype(first(points)))
  w = Vector{T}(undef, length(points))
  Threads.@threads for i in eachindex(points)
    w[i] = windfun(points[i])
  end
  w
end

winding(point::Point, mesh::Mesh) = winding((point,), mesh) |> first
