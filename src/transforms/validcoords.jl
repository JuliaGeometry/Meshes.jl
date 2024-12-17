# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ValidCoords(CRS)
    ValidCoords(code)

Retain the geometries within the domain of the
projection of type `CRS` or with EPSG/ESRI `code`.
"""
struct ValidCoords{CRS} <: CoordinateTransform end

ValidCoords(CRS) = ValidCoords{CRS}()

ValidCoords(code::Type{<:EPSG}) = ValidCoords(CoordRefSystems.get(code))

ValidCoords(code::Type{<:ESRI}) = ValidCoords(CoordRefSystems.get(code))

parameters(::ValidCoords{CRS}) where {CRS} = (; CRS)

preprocess(t::ValidCoords, d::Domain) = findall(g -> all(_indomain(t), pointify(g)), d)

function preprocess(t::ValidCoords, d::Mesh)
  indomain = _indomain(t)
  points = vertices(d)
  topo = topology(d)
  ∂₂₀ = Boundary{2,0}(topo)
  findall(1:nelements(d)) do elem
    is = ∂₂₀(elem)
    all(indomain(points[i]) for i in is)
  end
end

preprocess(t::ValidCoords, d::GeometrySet{<:Manifold,<:CRS,<:Union{Polytope,MultiPolytope}}) =
  findall(g -> all(_indomain(t), eachvertex(g)), d)

preprocess(t::ValidCoords, d::PointSet) = findall(_indomain(t), d)

apply(t::ValidCoords, d::Domain) = view(d, preprocess(t, d)), nothing

# -----------
# IO METHODS
# -----------

Base.show(io::IO, ::ValidCoords{CRS}) where {CRS} = print(io, "ValidCoords(CRS: $CRS)")

function Base.show(io::IO, ::MIME"text/plain", t::ValidCoords{CRS}) where {CRS}
  summary(io, t)
  println(io)
  print(io, "└─ CRS: $CRS")
end

# -----------------
# HELPER FUNCTIONS
# -----------------

_indomain(::ValidCoords{CRS}) where {CRS} = p -> indomain(CRS, coords(p))
