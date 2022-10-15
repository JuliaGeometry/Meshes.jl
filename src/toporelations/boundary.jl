# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Boundary{P,Q}(topology)

The boundary relation from rank `P` to smaller rank `Q` for
a given `topology`.
"""
struct Boundary{P,Q,D,T<:Topology} <: TopologicalRelation
  topology::T
end

function Boundary{P,Q}(topology) where {P,Q}
  D = paramdim(topology)
  T = typeof(topology)

  @assert D ≥ P > Q "invalid boundary relation"

  Boundary{P,Q,D,T}(topology)
end

# --------------
# GRID TOPOLOGY
# --------------

# quadrangle faces making up hexahedrons in 3D grid
function (∂::Boundary{3,2,3,T})(ind::Integer) where {T<:GridTopology}
  @error "not implemented"
end

# vertices of hexahedron on 3D grid
function (∂::Boundary{3,0,3,T})(ind::Integer) where {T<:GridTopology}
  t = ∂.topology
  cx, cy, cz = isperiodic(t)
  nx, ny, nz = size(t)

  i, j, k = elem2cart(t, ind)
  i₊ = cx ? mod1(i + 1, nx) : i + 1
  j₊ = cy ? mod1(j + 1, ny) : j + 1
  k₊ = cz ? mod1(k + 1, nz) : k + 1

  i1 = cart2corner(t, i , j , k )
  i2 = cart2corner(t, i₊, j , k )
  i3 = cart2corner(t, i₊, j₊, k )
  i4 = cart2corner(t, i , j₊, k )
  i5 = cart2corner(t, i , j , k₊)
  i6 = cart2corner(t, i₊, j , k₊)
  i7 = cart2corner(t, i₊, j₊, k₊)
  i8 = cart2corner(t, i , j₊, k₊)
  [i1, i2, i3, i4, i5, i6, i7, i8]
end

# vertices of quadrangle on 3D grid
function (∂::Boundary{2,0,3,T})(ind::Integer) where {T<:GridTopology}
  @error "not implemented"
end

# segment faces making up quadrangles in 2D grid
function (∂::Boundary{2,1,2,T})(ind::Integer) where {T<:GridTopology}
  t = ∂.topology
  cx, cy = isperiodic(t)
  nx, ny = size(t)

  # index in flipped grid
  i, j = elem2cart(t, ind)
  flip = GridTopology(ny, nx)
  find = cart2elem(flip, j, i)

  # vertical edges
  i1 = ind
  i2 = cy ? mod1(i1 + nx, nx*ny) : i1 + nx

  # horizontal edges
  i3 = find
  i4 = cx ? mod1(i3 + ny, nx*ny) : i3 + ny

  # horizontal edges come after vertical edges
  offset = nx*ny + !cy*nx
  i3 += offset
  i4 += offset

  [i1, i2, i3, i4]
end

# vertices of quadrangle on 2D grid
function (∂::Boundary{2,0,2,T})(ind::Integer) where {T<:GridTopology}
  t = ∂.topology
  cx, cy = isperiodic(t)
  nx, ny = size(t)

  i, j = elem2cart(t, ind)
  i₊ = cx ? mod1(i + 1, nx) : i + 1
  j₊ = cy ? mod1(j + 1, ny) : j + 1

  i1 = cart2corner(t, i , j )
  i2 = cart2corner(t, i₊, j )
  i3 = cart2corner(t, i₊, j₊)
  i4 = cart2corner(t, i , j₊)
  [i1, i2, i3, i4]
end

# vertices of segment on 2D grid
function (∂::Boundary{1,0,2,T})(ind::Integer) where {T<:GridTopology}
  t = ∂.topology
  cx, cy = isperiodic(t)
  nx, ny = size(t)

  @assert !(cx || cy) "not implemented"

  if ind ≤ nx*ny # vertical edges
    i1 = elem2corner(t, ind)
    i2 = i1 + 1
  elseif ind ≤ nx*(ny+1) # last vertical edges
    i, j = elem2cart(t, ind - nx)
    i1 = cart2corner(t, i, j+1)
    i2 = i1 + 1
  else # horizontal edges
    i1 = ind - nx*(ny+1)
    i, j = corner2cart(t, i1)
    i2 = cart2corner(t, i, j+1)
  end

  [i1, i2]
end

# vertices of segment on 1D grid
function (∂::Boundary{1,0,1,T})(ind::Integer) where {T<:GridTopology}
  t = ∂.topology
  c = first(isperiodic(t))
  n = first(size(t))

  i1 = ind
  i2 = c ? mod1(ind + 1, n) : ind + 1

  [i1, i2]
end

# -------------------
# HALF-EDGE TOPOLOGY
# -------------------

function (∂::Boundary{2,1,2,T})(elem::Integer) where {T<:HalfEdgeTopology}
  t = ∂.topology
  l = loop(half4elem(t, elem))
  v = CircularVector(l)
  [edge4pair(t, (v[i], v[i+1])) for i in 1:length(v)]
end

function (∂::Boundary{2,0,2,T})(elem::Integer) where {T<:HalfEdgeTopology}
  loop(half4elem(∂.topology, elem))
end

function (∂::Boundary{1,0,2,T})(edge::Integer) where {T<:HalfEdgeTopology}
  e = half4edge(∂.topology, edge)
  [e.head, e.half.head]
end

# ----------------
# SIMPLE TOPOLOGY
# ----------------

function (∂::Boundary{D,0,D,T})(ind::Integer) where {D,T<:SimpleTopology}
  collect(connec4elem(∂.topology, ind))
end