# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LengthUnits(unit)

Convert the length units of coordinates of a geometry or domain to `unit`.

## Examples

```julia
LengthUnits(u"cm")
LengthUnits(u"km")
```
"""
struct LengthUnits{U} <: CoordinateTransform
  unit::U
end

parameters(t::LengthUnits) = (; unit=t.unit)

applycoord(t::LengthUnits, v::Vec) = uconvert.(t.unit, v)

function applycoord(t::LengthUnits, p::Point{<:Any,<:Cartesian})
  c = CoordRefSystems.cvalues(coords(p))
  Point(Cartesian{datum(crs(p))}(uconvert.(t.unit, c)))
end

function applycoord(t::LengthUnits, p::Point{<:Any,<:Polar})
  c = coords(p)
  ρ = uconvert(t.unit, c.ρ)
  Point(Polar{datum(crs(p))}(ρ, c.ϕ))
end

function applycoord(t::LengthUnits, p::Point{<:Any,<:Cylindrical})
  c = coords(p)
  ρ = uconvert(t.unit, c.ρ)
  z = uconvert(t.unit, c.z)
  Point(Cylindrical{datum(crs(p))}(ρ, c.ϕ, z))
end

function applycoord(t::LengthUnits, p::Point{<:Any,<:Spherical})
  c = coords(p)
  r = uconvert(t.unit, c.r)
  Point(Spherical{datum(crs(p))}(r, c.θ, c.ϕ))
end

applycoord(::LengthUnits, ::Point) = throw(ArgumentError("only the length units of Basic CRSs can be changed"))

# --------------
# SPECIAL CASES
# --------------

applycoord(t::LengthUnits, g::RectilinearGrid) = RectilinearGrid{datum(crs(g))}(map(x -> uconvert.(t.unit, x), xyz(g)))

applycoord(t::LengthUnits, g::StructuredGrid) = StructuredGrid{datum(crs(g))}(map(X -> uconvert.(t.unit, X), XYZ(g)))
