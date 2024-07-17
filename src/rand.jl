# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    rand([rng], G, crs=Cartesian3D)

Generate a random geometry of type `G` with CRS `crs`,
optionally passing a random number generator `rng`.

# Examples

```julia
rand(Point)
rand(Triangle)
rand(Point, crs=Cartesian2D)
rand(Triangle, crs=LatLon)
```
"""
Random.rand(G::Type{<:Geometry}; kwargs...) = rand(Random.default_rng(), G; kwargs...)
Random.rand(rng::Random.AbstractRNG, G::Type{<:Geometry}; crs=Cartesian3D) = _rand(rng, G, crs)

"""
    rand([rng], G, n, crs=Cartesian3D)

Generate a vector of `n` random geometries of type `G` with CRS `crs`,
optionally passing a random number generator `rng`.

# Examples

```julia
rand(Point, 10)
rand(Triangle, 10)
rand(Point, 10, crs=Cartesian2D)
rand(Triangle, 10, crs=LatLon)
```
"""
Random.rand(G::Type{<:Geometry}, n::Int; kwargs...) = rand(Random.default_rng(), G, n; kwargs...)
Random.rand(rng::Random.AbstractRNG, G::Type{<:Geometry}, n::Int; kwargs...) = [rand(rng, G; kwargs...) for _ in 1:n]

# ----------------
# IMPLEMENTATIONS
# ----------------

_rand(rng::Random.AbstractRNG, ::Type{Point}, CRS) = Point(rand(rng, CRS))

_rand(rng::Random.AbstractRNG, ::Type{Ray}, CRS) = Ray(_rand(rng, Point, CRS), _rvec(rng, CRS))

_rand(rng::Random.AbstractRNG, ::Type{Line}, CRS) = Line(_rand(rng, Point, CRS), _rand(rng, Point, CRS))

_rand(rng::Random.AbstractRNG, ::Type{BezierCurve}, CRS) = BezierCurve(_rvector(rng, CRS, 5))

_rand(rng::Random.AbstractRNG, ::Type{Plane}, CRS) = Plane(_rand(rng, Point, CRS), _rvec(rng, CRS))

function _rand(rng::Random.AbstractRNG, ::Type{Box}, CRS)
  min = _rand(rng, Point, CRS)
  max = min + _rvec(rng, CRS)
  Box(min, max)
end

_rand(rng::Random.AbstractRNG, ::Type{Ball}, CRS) = Ball(_rand(rng, Point, CRS), _rlen(rng))

_rand(rng::Random.AbstractRNG, ::Type{Sphere}, CRS) = Sphere(_rand(rng, Point, CRS), _rlen(rng))

_rand(rng::Random.AbstractRNG, ::Type{Ellipsoid}, CRS) =
  Ellipsoid((_rlen(rng), _rlen(rng), _rlen(rng)), _rand(rng, Point, CRS), rand(rng, QuatRotation))

_rand(rng::Random.AbstractRNG, ::Type{Disk}, CRS) = Disk(_rand(rng, Plane, CRS), _rlen(rng))

_rand(rng::Random.AbstractRNG, ::Type{Circle}, CRS) = Circle(_rand(rng, Plane, CRS), _rlen(rng))

_rand(rng::Random.AbstractRNG, ::Type{Cylinder}, CRS) =
  Cylinder(_rand(rng, Plane, CRS), _rand(rng, Plane, CRS), _rlen(rng))

_rand(rng::Random.AbstractRNG, ::Type{CylinderSurface}, CRS) =
  CylinderSurface(_rand(rng, Plane, CRS), _rand(rng, Plane, CRS), _rlen(rng))

_rand(rng::Random.AbstractRNG, ::Type{Cone}, CRS) = Cone(_rand(rng, Disk, CRS), _rand(rng, Point, CRS))

_rand(rng::Random.AbstractRNG, ::Type{ConeSurface}, CRS) = ConeSurface(_rand(rng, Disk, CRS), _rand(rng, Point, CRS))

function _rand(rng::Random.AbstractRNG, ::Type{Frustum}, CRS)
  bottom = _rand(rng, Disk, CRS)
  ax = normal(plane(bottom))
  topplane = Plane(center(bottom) + rand(rng) * ax, ax)
  top = Disk(topplane, _rlen(rng))
  Frustum(bottom, top)
end

function _rand(rng::Random.AbstractRNG, ::Type{FrustumSurface}, CRS)
  bottom = _rand(rng, Disk, CRS)
  ax = normal(plane(bottom))
  topplane = Plane(center(bottom) + rand(rng) * ax, ax)
  top = Disk(topplane, _rlen(rng))
  FrustumSurface(bottom, top)
end

_rand(rng::Random.AbstractRNG, ::Type{ParaboloidSurface}, CRS) =
  ParaboloidSurface(_rand(rng, Point, CRS), _rlen(rng), _rlen(rng))

_rand(rng::Random.AbstractRNG, ::Type{Torus}, CRS) =
  Torus(_rand(rng, Point, CRS), _rvec(rng, CRS), _rlen(rng), _rlen(rng))

_rand(rng::Random.AbstractRNG, ::Type{Segment}, CRS) = Segment(_rtuple(rng, CRS, 2))

_rand(rng::Random.AbstractRNG, ::Type{Rope}, CRS) = Rope(_rvector(rng, CRS, rand(rng, 2:50)))

function _rand(rng::Random.AbstractRNG, ::Type{Ring}, CRS)
  v = _rvector(rng, CRS, rand(rng, 3:50))
  while first(v) == last(v)
    v = _rvector(rng, CRS, rand(rng, 3:50))
  end
  Ring(v)
end

_rand(rng::Random.AbstractRNG, ::Type{Ngon{N}}, CRS) where {N} = Ngon{N}(_rtuple(rng, CRS, N))

_rand(rng::Random.AbstractRNG, ::Type{PolyArea}, CRS) = PolyArea(_rand(rng, Ring, CRS))

_rand(rng::Random.AbstractRNG, ::Type{Tetrahedron}, CRS) = Tetrahedron(_rtuple(rng, CRS, 4))

_rand(rng::Random.AbstractRNG, ::Type{Hexahedron}, CRS) = Hexahedron(_rtuple(rng, CRS, 8))

_rand(rng::Random.AbstractRNG, ::Type{Pyramid}, CRS) = Pyramid(_rtuple(rng, CRS, 5))

_rand(rng::Random.AbstractRNG, ::Type{Wedge}, CRS) = Wedge(_rtuple(rng, CRS, 6))

# -----------------
# HELPER FUNCTIONS
# -----------------

_rvector(rng::Random.AbstractRNG, CRS, n) = [_rand(rng, Point, CRS) for _ in 1:n]

_rtuple(rng::Random.AbstractRNG, CRS, n) = ntuple(_ -> _rand(rng, Point, CRS), n)

_rvec(rng::Random.AbstractRNG, CRS) = rand(rng, Vec{CoordRefSystems.ndims(CRS),Met{Float64}})

_rlen(rng::Random.AbstractRNG) = rand(rng, Met{Float64})
