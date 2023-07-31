# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

Base.issubset(p::Point, g::Geometry) = p ∈ g

Base.issubset(s₁::Segment, s₂::Segment) = all(∈(s₂), vertices(s₁))

Base.issubset(s::Segment, b::Box) = all(∈(b), vertices(s))

Base.issubset(s::Segment, b::Ball) = all(∈(b), vertices(s))

Base.issubset(t₁::Triangle, t₂::Triangle) = all(∈(t₂), vertices(t₁))

Base.issubset(t::Triangle, p::Polygon) =
  all(∈(p), vertices(t)) && (!hasintersect(boundary(t), boundary(p)) || vertices(t) ⊆ boundary(p))

Base.issubset(t::Triangle, b::Box) = all(∈(b), vertices(t))

Base.issubset(t::Triangle, b::Ball) = all(∈(b), vertices(t))

Base.issubset(t₁::Tetrahedron, t₂::Tetrahedron) = all(∈(t₂), vertices(t₁))

Base.issubset(t::Tetrahedron, b::Box) = all(∈(b), vertices(t))

Base.issubset(t::Tetrahedron, b::Ball) = all(∈(b), vertices(t))

Base.issubset(g₁::Geometry, g₂::Geometry) = all(g -> g ⊆ g₂, simplexify(g₁))
