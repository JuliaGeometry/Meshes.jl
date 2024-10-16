# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Proj(CRS; boundary=false)
    Proj(code; boundary=false)

Convert the coordinates of geometry or domain to a given
coordinate reference system `CRS` or EPSG/ESRI `code`.

Optionally, the transform samples the `boundary` of
polytopes, if this option is `true`, to handle distortions
that occur in manifold conversions.

## Examples

```julia
Proj(Polar)
Proj(WebMercator)
Proj(Mercator{WGS84Latest})
Proj(EPSG{3395})
Proj(ESRI{54017})
Proj(Robinson, boundary=true)
```
"""
struct Proj{CRS,Boundary} <: CoordinateTransform end

Proj(CRS; boundary=false) = Proj{CRS,boundary}()

Proj(code::Type{<:EPSG}; kwargs...) = Proj(CoordRefSystems.get(code); kwargs...)

Proj(code::Type{<:ESRI}; kwargs...) = Proj(CoordRefSystems.get(code); kwargs...)

parameters(::Proj{CRS,Boundary}) where {CRS,Boundary} = (CRS=CRS, boundary=Boundary)

# avoid constructing a new geometry or domain when the CRS is the same
function apply(t::Proj{CRS}, g::GeometryOrDomain) where {CRS}
  g′ = crs(g) <: CRS ? g : applycoord(t, g)
  g′, nothing
end

# convert the CRS and preserve the manifold
applycoord(::Proj{CRS}, p::Point{<:🌐}) where {CRS<:Basic} = Point{🌐}(convert(CRS, coords(p)))

# convert the CRS and (possibly) change the manifold
applycoord(::Proj{CRS}, p::Point{<:🌐}) where {CRS<:Projected} = _proj(CRS, p)
applycoord(::Proj{CRS}, p::Point{<:🌐}) where {CRS<:Geographic} = _proj(CRS, p)
applycoord(::Proj{CRS}, p::Point{<:𝔼}) where {CRS<:Basic} = _proj(CRS, p)
applycoord(::Proj{CRS}, p::Point{<:𝔼}) where {CRS<:Projected} = _proj(CRS, p)
applycoord(::Proj{CRS}, p::Point{<:𝔼}) where {CRS<:Geographic} = _proj(CRS, p)

applycoord(::Proj, v::Vec) = v

# --------------
# SPECIAL CASES
# --------------

applycoord(t::Proj{<:Projected}, g::Primitive{<:🌐}) = TransformedGeometry(g, t)

applycoord(t::Proj{<:Geographic}, g::Primitive{<:𝔼}) = TransformedGeometry(g, t)

applycoord(t::Proj{<:Projected,true}, g::Polytope{K,<:🌐}) where {K} = TransformedGeometry(g, t)

applycoord(t::Proj{<:Geographic,true}, g::Polytope{K,<:𝔼}) where {K} = TransformedGeometry(g, t)

applycoord(t::Proj, g::RegularGrid) = TransformedGrid(g, t)

applycoord(t::Proj, g::RectilinearGrid) = TransformedGrid(g, t)

applycoord(t::Proj, g::StructuredGrid) = TransformedGrid(g, t)

# -----------
# IO METHODS
# -----------

Base.show(io::IO, ::Proj{CRS,Boundary}) where {CRS,Boundary} = print(io, "Proj(CRS: $CRS, boundary: $Boundary)")

function Base.show(io::IO, ::MIME"text/plain", t::Proj{CRS,Boundary}) where {CRS,Boundary}
  summary(io, t)
  println(io)
  println(io, "├─ CRS: $CRS")
  print(io, "└─ boundary: $Boundary")
end

# -----------------
# HELPER FUNCTIONS
# -----------------

_proj(CRS, p) = Point(convert(CRS, coords(p)))
