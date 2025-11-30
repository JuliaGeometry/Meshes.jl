# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    geometry₁ ⊆ geometry₂

Tells whether or not `geometry₁` is contained in `geometry₂`.
"""
function Base.issubset(g₁::Geometry, g₂::Geometry)
  if g₁ == g₂
    return true
  elseif isconvex(g₁) && isconvex(g₂)
    return _boundarypoints(g₁) ⊆ g₂
  elseif isconvex(g₁)
    return _boundarypoints(g₁) ⊆ g₂ && !intersects(g₁, !g₂)
  else
    return all(g -> g ⊆ g₂, simplexify(g₁))
  end
end

_boundarypoints(p::Primitive) = boundarypoints(p)

_boundarypoints(p::Polytope) = vertices(p)

# --------------
# OPTIMIZATIONS
# --------------

Base.issubset(p::Point, g::Geometry) = p ∈ g

Base.issubset(b₁::Box, b₂::Box) = minimum(b₁) ∈ b₂ && maximum(b₁) ∈ b₂

# ---------------
# IMPLEMENTATION
# ---------------

Base.issubset(points, geom::Geometry) = all(∈(geom), points)
