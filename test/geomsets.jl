@testitem "GeometrySet" setup = [Setup] begin
  s = Segment(cart(0, 0), cart(1, 1))
  t = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
  p = PolyArea(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
  gset = GeometrySet([s, t, p])
  @test crs(gset) <: Cartesian{NoDatum}
  @test Meshes.lentype(gset) == ℳ
  @test [centroid(gset, i) for i in 1:3] == cart.([(1 / 2, 1 / 2), (1 / 3, 1 / 3), (1 / 2, 1 / 2)])

  s = Segment(cart(0, 0), cart(1, 1))
  t = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
  geoms = [s, t]
  gset1 = GeometrySet(geoms)
  gset2 = GeometrySet(g for g in geoms)
  @test gset1 == gset2
  @test parent(gset1) === geoms

  # make sure that eltype is inferred properly
  # https://github.com/JuliaGeometry/Meshes.jl/issues/643
  geoms = Vector{Segment}()
  push!(geoms, Segment(cart(0, 0), cart(1, 0)))
  push!(geoms, Segment(cart(1, 0), cart(1, 1)))
  push!(geoms, Segment(cart(1, 1), cart(0, 0)))
  gset = GeometrySet(geoms)
  @test eltype(gset) <: Segment
  @test nelements(gset) == 3

  # construction with custom vectors of geometries
  # https://github.com/JuliaEarth/GeoStats.jl/issues/551
  geoms = @SVector [cart(0, 0), cart(1, 1)]
  gset = GeometrySet(geoms)
  @test eltype(gset) <: Point
  @test nelements(gset) == 2

  # different CRS
  s = Segment(latlon(0, 0), latlon(1, 1))
  t = Triangle(latlon(0, 0), latlon(0, 1), latlon(1, 0)) |> Proj(PlateCarree)
  q = Quadrangle(latlon(0, 0), latlon(0, 1), latlon(1, 1), latlon(1, 0)) |> Proj(WebMercator)
  geoms = [s, t, q]
  gset = GeometrySet(geoms)
  @test eltype(gset) <: Polytope
  @test crs(gset) <: LatLon
  gset = GeometrySet(g for g in geoms)
  @test eltype(gset) <: Polytope
  @test crs(gset) <: LatLon
  geoms = [t, s, q]
  gset = GeometrySet(geoms)
  @test eltype(gset) <: Polytope
  @test crs(gset) <: PlateCarree
  gset = GeometrySet(g for g in geoms)
  @test eltype(gset) <: Polytope
  @test crs(gset) <: PlateCarree
  geoms = [q, s, t]
  gset = GeometrySet(geoms)
  @test eltype(gset) <: Polytope
  @test crs(gset) <: WebMercator
  gset = GeometrySet(g for g in geoms)
  @test eltype(gset) <: Polytope
  @test crs(gset) <: WebMercator

  # conversion
  grid = cartgrid(10, 10)
  gset = convert(GeometrySet, grid)
  @test gset isa GeometrySet
  @test nelements(gset) == 100
  @test eltype(gset) <: Quadrangle
end

@testitem "PointSet" setup = [Setup] begin
  pset = PointSet([cart(1), cart(2), cart(3)])
  @test embeddim(pset) == 1
  @test crs(pset) <: Cartesian{NoDatum}
  @test Meshes.lentype(pset) === ℳ
  @test nelements(pset) == 3
  @test eltype(pset) <: Point

  pset = PointSet([cart(1, 1), cart(2, 2), cart(3, 3)])
  @test embeddim(pset) == 2
  @test crs(pset) <: Cartesian{NoDatum}
  @test Meshes.lentype(pset) === ℳ
  @test nelements(pset) == 3
  @test eltype(pset) <: Point

  pset = PointSet([cart(1, 1, 1), cart(2, 2, 2), cart(3, 3, 3)])
  @test embeddim(pset) == 3
  @test crs(pset) <: Cartesian{NoDatum}
  @test Meshes.lentype(pset) === ℳ
  @test nelements(pset) == 3
  @test eltype(pset) <: Point

  pset1 = PointSet([cart(1, 2, 3), cart(4, 5, 6)])
  pset2 = PointSet(cart(1, 2, 3), cart(4, 5, 6))
  pset3 = PointSet([T.((1, 2, 3)), T.((4, 5, 6))])
  pset4 = PointSet(T.((1, 2, 3)), T.((4, 5, 6)))
  @test pset1 == pset2 == pset3 == pset4
  for pset in [pset1, pset2, pset3, pset4]
    @test embeddim(pset) == 3
    @test Meshes.lentype(pset) === ℳ
    @test nelements(pset) == 2
    @test pset[1] == cart(1, 2, 3)
    @test pset[2] == cart(4, 5, 6)
  end

  pset = PointSet(cart.([(0, 0), (1, 0), (0, 1)]))
  @test centroid(pset) == cart(1 / 3, 1 / 3)

  pset = PointSet(cart.([(1, 0), (0, 1)]))
  @test nelements(pset) == 2
  @test centroid(pset, 1) == cart(1, 0)
  @test centroid(pset, 2) == cart(0, 1)

  pset = PointSet(cart.([(0, 0), (1, 0), (0, 1)]))
  @test measure(pset) == zero(T) * u"m"

  # constructor with iterator
  points = cart.([(1, 0), (0, 1)])
  pset1 = PointSet(points)
  pset2 = PointSet(p for p in points)
  @test pset1 == pset2

  # CRS propagation
  pset = PointSet(merc.([(0, 0), (1, 0), (0, 1)]))
  @test crs(centroid(pset)) === crs(pset)

  pset = PointSet(cart.([(1, 0), (0, 1)]))
  @test sprint(show, pset) == "2 PointSet"
  if T === Float32
    @test sprint(show, MIME"text/plain"(), pset) == """
    2 PointSet
    ├─ Point(x: 1.0f0 m, y: 0.0f0 m)
    └─ Point(x: 0.0f0 m, y: 1.0f0 m)"""
  elseif T == Float64
    @test sprint(show, MIME"text/plain"(), pset) == """
    2 PointSet
    ├─ Point(x: 1.0 m, y: 0.0 m)
    └─ Point(x: 0.0 m, y: 1.0 m)"""
  end
end
