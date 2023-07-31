# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    isconvex(geometry)

Tells whether or not the `geometry` is convex.
"""
function isconvex end

isconvex(g::Geometry) = isconvex(typeof(g))

isconvex(::Type{<:Point}) = true

isconvex(::Type{<:Segment}) = true

isconvex(::Type{<:Ray}) = true

isconvex(::Type{<:Line}) = true

isconvex(::Type{<:Plane}) = true

isconvex(::Type{<:Box}) = true

isconvex(::Type{<:Ball}) = true

isconvex(::Type{<:Sphere}) = false

isconvex(::Type{<:Disk}) = true

isconvex(::Type{<:Circle}) = false

isconvex(::Type{<:Cone}) = true

isconvex(::Type{<:ConeSurface}) = false

isconvex(::Type{<:Cylinder}) = true

isconvex(::Type{<:CylinderSurface}) = false

isconvex(::Type{<:Torus}) = false

isconvex(p::Polygon{Dim,T}) where {Dim,T} = issimple(p) && all(≤(T(π)), innerangles(boundary(p)))

isconvex(::Type{<:Triangle}) = true

isconvex(::Type{<:Tetrahedron}) = true

# --------------
# OPTIMIZATIONS
# --------------

isconvex(::Triangle) = true

function isconvex(q::Quadrangle{2})
  v = vertices(q)
  d1 = Segment(v[1], v[3])
  d2 = Segment(v[2], v[4])
  intersects(d1, d2)
end

function isconvex(q::Quadrangle{3})
  v = vertices(q)
  iscoplanar(v...) &&
  isconvex(Quadrangle(proj2D(collect(v))...))
end
