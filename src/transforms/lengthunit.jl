# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LengthUnit(unit)

Convert the length unit of coordinates of a geometry or domain to `unit`.

## Examples

```julia
LengthUnit(u"cm")
LengthUnit(u"km")
```
"""
struct LengthUnit{U} <: CoordinateTransform
  unit::U
end

parameters(t::LengthUnit) = (; unit=t.unit)

applycoord(t::LengthUnit, p::Point) = Point(_lenunit(coords(p), t.unit))

applycoord(t::LengthUnit, v::Vec) = uconvert.(t.unit, v)

# --------------
# SPECIAL CASES
# --------------

applycoord(t::LengthUnit, len::Len) = uconvert(t.unit, len)

applycoord(t::LengthUnit, lens::NTuple{Dim,Len}) where {Dim} = uconvert.(t.unit, lens)

function applycoord(t::LengthUnit, g::RegularGrid)
  dims = size(g)
  orig = applycoord(t, minimum(g))
  spac = map(s -> applycoord(t, s), spacing(g))
  offs = offset(g)
  RegularGrid(dims, orig, spac, offs)
end

applycoord(t::LengthUnit, g::RectilinearGrid) = TransformedGrid(g, t)

applycoord(t::LengthUnit, g::StructuredGrid) = TransformedGrid(g, t)

# -----------------
# HELPER FUNCTIONS
# -----------------

function _lenunit(c::Cartesian, u)
  d = datum(c)
  v = CoordRefSystems.values(c)
  Cartesian{d}(uconvert.(u, v))
end

function _lenunit(c::Polar, u)
  d = datum(c)
  ρ = uconvert(u, c.ρ)
  Polar{d}(ρ, c.ϕ)
end

function _lenunit(c::Cylindrical, u)
  d = datum(c)
  ρ = uconvert(u, c.ρ)
  z = uconvert(u, c.z)
  Cylindrical{d}(ρ, c.ϕ, z)
end

function _lenunit(c::Spherical, u)
  d = datum(c)
  r = uconvert(u, c.r)
  Spherical{d}(r, c.θ, c.ϕ)
end

_lenunit(c, _) = throw(ArgumentError("the length unit of $(prettyname(c)) cannot be changed"))
