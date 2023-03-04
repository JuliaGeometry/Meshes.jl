# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    merge(object₁, object₂)

Merge `object₁` with `object₂`, i.e. concatenate
the vertices and adjust the connectivities accordingly.
"""
function Base.merge(m₁::Mesh, m₂::Mesh)
  v₁ = vertices(m₁)
  v₂ = vertices(m₂)
  t₁ = topology(m₁)
  t₂ = topology(m₂)

  # concatenate vertices
  points = [v₁; v₂]

  # concatenate connectivities
  offset = length(v₁)
  connec₁ = collect(elements(t₁))
  connec₂ = map(elements(t₂)) do e
    PL = pltype(e)
    c = indices(e)
    c′ = ntuple(i -> c[i] + offset, length(c))
    connect(c′, PL)
  end
  connec = [connec₁; connec₂]

  return SimpleMesh(points, connec)
end

Base.merge(g::Geometry, m::Mesh) = merge(discretize(g), m)

Base.merge(m::Mesh, g::Geometry) = merge(m, discretize(g))

Base.merge(g₁::Geometry, g₂::Geometry) = Multi([g₁, g₂])

Base.merge(m₁::Multi, m₂::Multi) = Multi([collect(m₁); collect(m₂)])

Base.merge(m::Multi, g::Geometry) = Multi([collect(m); g])

Base.merge(g::Geometry, m::Multi) = Multi([g; collect(m)])
