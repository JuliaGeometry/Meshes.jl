# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    InDomain(CRS)
    InDomain(code)

Retain the geometries within the domain of the
projection of type `CRS` or with EPSG/ESRI `code`.
"""
struct InDomain{CRS} <: CoordinateTransform end

InDomain(CRS) = InDomain{CRS}()

InDomain(code::Type{<:EPSG}) = InDomain(CoordRefSystems.get(code))

InDomain(code::Type{<:ESRI}) = InDomain(CoordRefSystems.get(code))

parameters(::InDomain{CRS}) where {CRS} = (; CRS)

preprocess(t::InDomain, d::Domain) = findall(g -> all(_indomain(t), pointify(g)), d)

function preprocess(t::InDomain, d::Mesh)
  indomain = _indomain(t)
  points = vertices(d)
  topo = topology(d)
  ∂₂₀ = Boundary{2,0}(topo)
  findall(1:nelements(d)) do elem
    is = ∂₂₀(elem)
    all(indomain(points[i]) for i in is)
  end
end

preprocess(t::InDomain, d::GeometrySet{<:Manifold,<:CRS,<:Union{Polytope,MultiPolytope}}) =
  findall(g -> all(_indomain(t), eachvertex(g)), d)

preprocess(t::InDomain, d::PointSet) = findall(_indomain(t), d)

apply(t::InDomain, d::Domain) = view(d, preprocess(t, d)), nothing

# -----------
# IO METHODS
# -----------

Base.show(io::IO, ::InDomain{CRS}) where {CRS} = print(io, "InDomain(CRS: $CRS)")

function Base.show(io::IO, ::MIME"text/plain", t::InDomain{CRS}) where {CRS}
  summary(io, t)
  println(io)
  print(io, "└─ CRS: $CRS")
end

# -----------------
# HELPER FUNCTIONS
# -----------------

_indomain(::InDomain{CRS}) where {CRS} = p -> indomain(CRS, coords(p))
