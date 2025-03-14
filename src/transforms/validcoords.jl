# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ValidCoords(CRS)

Retain the geometries within the domain of the projected `CRS`.

    ValidCoords(code)

Alternatively, specify the projected `CRS` using a EPSG/ESRI `code`.
"""
struct ValidCoords{CRS} <: GeometricTransform end

ValidCoords(CRS) = ValidCoords{CRS}()

ValidCoords(code::Type{<:EPSG}) = ValidCoords(CoordRefSystems.get(code))

ValidCoords(code::Type{<:ESRI}) = ValidCoords(CoordRefSystems.get(code))

parameters(::ValidCoords{CRS}) where {CRS} = (; CRS)

preprocess(t::ValidCoords, d::Domain) = findall(g -> _isvalid(t, g), d)

function preprocess(t::ValidCoords, d::Mesh)
  points = vertices(d)
  topo = topology(d)
  findall(elements(topo)) do elem
    is = indices(elem)
    all(_isvalid(t, points[i]) for i in is)
  end
end

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

_isvalid(::ValidCoords{CRS}, p::Point) where {CRS} = indomain(CRS, coords(p))
_isvalid(t::ValidCoords, g::Polytope) = all(p -> _isvalid(t, p), eachvertex(g))
_isvalid(t::ValidCoords, g::MultiPolytope) = all(p -> _isvalid(t, p), eachvertex(g))
_isvalid(t::ValidCoords, g::Geometry) = all(p -> _isvalid(t, p), pointify(g))
